"""True-orientation PDF text geometry; no drawing and no relaxed thresholds."""
from __future__ import annotations
import math
from itertools import combinations


def area(poly):
    return abs(sum(a[0]*b[1]-b[0]*a[1] for a,b in zip(poly,poly[1:]+poly[:1])))/2


def cross(a,b,c):
    return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])


def intersection(subject,clip):
    sign=1 if sum(a[0]*b[1]-b[0]*a[1] for a,b in zip(clip,clip[1:]+clip[:1]))>=0 else -1
    out=list(subject)
    for a,b in zip(clip,clip[1:]+clip[:1]):
        inp=out;out=[]
        if not inp:break
        previous=inp[-1];dp=sign*cross(a,b,previous)
        for current in inp:
            dc=sign*cross(a,b,current)
            if (dc>=-1e-9)!=(dp>=-1e-9):
                t=dp/(dp-dc)
                out.append((previous[0]+t*(current[0]-previous[0]),previous[1]+t*(current[1]-previous[1])))
            if dc>=-1e-9:out.append(current)
            previous=current;dp=dc
    return out


def quad_from_chars(pymupdf,direction,span,chars):
    u=tuple(direction);v=(-u[1],u[0]);points=[]
    for c in chars:
        q=pymupdf.recover_char_quad(direction,span,c)
        points.extend([(p.x,p.y) for p in (q.ul,q.ur,q.lr,q.ll)])
    us=[x*u[0]+y*u[1] for x,y in points];vs=[x*v[0]+y*v[1] for x,y in points]
    bounds=(min(us),min(vs),max(us),max(vs))
    lo,bot,hi,top=bounds
    poly=[(a*u[0]+b*v[0],a*u[1]+b*v[1]) for a,b in ((lo,bot),(hi,bot),(hi,top),(lo,top))]
    return dict(poly=poly,u=u,v=v,bounds=bounds)


def install_oriented_audit(mod):
    original=mod.audit_page_geometry

    def audit(page,**kwargs):
        old,info=original(page,**kwargs)
        # All non-text and clipping categories retain the original engine.
        findings=[f for f in old if f.kind not in ('text-text','text-stroke')]
        inset=kwargs.get('text_inset_pt',.6);threshold=kwargs.get('min_overlap_ratio',.05)
        shapes=page.oriented_text
        for a,b in combinations(page.texts,2):
            if mod.rect_intersection(a.bbox,b.bbox) is None:continue
            pa=shapes[a.index]['poly'];pb=shapes[b.index]['poly']
            ratio=area(intersection(pa,pb))/min(area(pa),area(pb))
            if ratio>=threshold:
                findings.append(mod.CollisionFinding(severity='FAIL',kind='text-text',page=page.page,
                    message=f'Oriented text quadrilaterals overlap by {ratio:.1%} of the smaller label',
                    text=a.text,text_bbox=a.bbox,other_text=b.text,other_bbox=b.bbox))
        for t in page.texts:
            shape=shapes[t.index];u=shape['u'];v=shape['v'];lo,bot,hi,top=shape['bounds']
            trim=min(max(0,inset),.2*min(hi-lo,top-bot));hits=[];rects=[]
            for stroke in page.strokes:
                if mod.rect_intersection(t.bbox,mod.rect_expand(stroke.bbox,stroke.width/2)) is None:continue
                pad=max(0,stroke.width)/2
                local=(lo+trim-pad,bot+trim-pad,hi-trim+pad,top-trim+pad)
                def project(p):return (p[0]*u[0]+p[1]*u[1],p[0]*v[0]+p[1]*v[1])
                if any(mod.segment_intersects_rect((project(a),project(b)),local) for a,b in stroke.segments):
                    hits.append(stroke.index);rects.append(stroke.bbox)
            if hits:
                findings.append(mod.CollisionFinding(severity='FAIL',kind='text-stroke',page=page.page,
                    message=f'Oriented text interior crossed by {len(hits)} stroked paths',text=t.text,
                    text_bbox=t.bbox,other_bbox=mod.union_rects(rects),object_count=len(hits),object_indexes=tuple(hits)))
        return sorted(findings,key=mod._finding_sort_key),info

    mod.audit_page_geometry=audit


def self_test():
    def rectangle(cx,cy,length=12,height=1,angle=math.pi/4):
        u=(math.cos(angle),math.sin(angle));v=(-u[1],u[0])
        return [(cx+x*u[0]+y*v[0],cy+x*u[1]+y*v[1]) for x,y in
                ((-length/2,-height/2),(length/2,-height/2),(length/2,height/2),(-length/2,height/2))]
    a=rectangle(0,0);b=rectangle(-1.5/math.sqrt(2),1.5/math.sqrt(2))
    # Parallel rotated labels have heavily intersecting horizontal bboxes but no true overlap.
    assert area(intersection(a,b))<1e-9
    c=rectangle(-.5/math.sqrt(2),.5/math.sqrt(2))
    assert abs(area(intersection(a,c))/area(a)-.5)<1e-9
    assert abs(area(intersection(a,a))-12)<1e-9
    assert area(intersection(a,rectangle(50,50)))==0
    assert abs(area(intersection([(0,0),(2,0),(2,1),(0,1)],[(1,0),(3,0),(3,1),(1,1)]))-1)<1e-9
    return 'PASS: disjoint rotated labels, true rotated overlap, identity and axis-aligned calibration'


if __name__=='__main__':print(self_test())
