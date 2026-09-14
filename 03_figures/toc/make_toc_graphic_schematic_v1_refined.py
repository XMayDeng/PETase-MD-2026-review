#!/usr/bin/env python3
"""Author-requested V1 arrow cleanup, group spacing and card centering.

Remove the central L10-to-L20 arrow and expose the shafts of the two grey
reading-order arrows. Preserve shapes, sizes and colors while moving
the L10/L20 groups closer and centering content inside the two cards.
Use the author-approved heading "Interaction geometry" for the right panel.
Align the three graphical groups and both arrows on a shared vertical center.
Use matching typography and vertical placement for the two panel headings.
Render from native drawing objects, never retouch raster pixels.
"""
from pathlib import Path
import hashlib
import json

import numpy as np
from PIL import Image
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch
from matplotlib.transforms import Affine2D, Bbox
import make_toc_graphic_schematic as original

HERE=Path(__file__).resolve().parent
OUTPUT=HERE.parents[1]/'figures'/'TOC_graphic.png'


def center_card_contents(figure):
    ax=figure.axes[0]
    figure.canvas.draw()
    renderer=figure.canvas.get_renderer()
    boxes=[p for p in ax.patches if isinstance(p,FancyBboxPatch)]
    assert len(boxes)==2
    shifts=[]
    for box in boxes:
        y0=box.get_y();y1=y0+box.get_height()
        content=[t for t in ax.texts if t.get_position()[0]>=.72 and y0<t.get_position()[1]<y1]
        content.extend(line for line in ax.lines if min(line.get_xdata())>=.72 and y0<np.mean(line.get_ydata())<y1)
        for patch in ax.patches:
            if isinstance(patch,FancyBboxPatch):continue
            if isinstance(patch,original.Circle):x,y=patch.center
            elif isinstance(patch,original.Polygon):x,y=np.mean(patch.get_xy(),axis=0)
            else:continue
            if x>=.72 and y0<y<y1:content.append(patch)
        assert len(content)==(10 if y0>.4 else 6),(y0,len(content))
        bounds=Bbox.union([a.get_window_extent(renderer) for a in content])
        target=box.get_window_extent(renderer)
        delta=np.array([(target.x0+target.x1-bounds.x0-bounds.x1)/2,
                        (target.y0+target.y1-bounds.y0-bounds.y1)/2])
        shift=Affine2D().translate(*delta)
        for artist in content:
            transform=artist.get_data_transform() if isinstance(artist,Patch) else artist.get_transform()
            artist.set_transform(transform+shift)
        figure.canvas.draw()
        after=Bbox.union([a.get_window_extent(figure.canvas.get_renderer()) for a in content])
        assert np.allclose(after.get_points(),bounds.get_points()+delta)
        shifts.append({'card':'TfCut1' if y0>.4 else 'HiC',
                       'dx_pixels':float(delta[0]),'dy_pixels_down':float(-delta[1])})
    return shifts


def tighten_center_groups(figure):
    ax=figure.axes[0]
    figure.canvas.draw()
    renderer=figure.canvas.get_renderer()
    groups={'L10':[],'L20':[]}
    for artist in [*ax.patches,*ax.lines]:
        if isinstance(artist,(FancyArrowPatch,FancyBboxPatch)):continue
        bounds=artist.get_window_extent(renderer)
        x=ax.transData.inverted().transform(((bounds.x0+bounds.x1)/2,
                                             (bounds.y0+bounds.y1)/2))[0]
        if .30<x<.69:groups['L10' if x<.48 else 'L20'].append(artist)
    for text in ax.texts:
        if text.get_text() in groups:groups[text.get_text()].append(text)
    assert {key:len(value) for key,value in groups.items()}=={'L10':10,'L20':18}
    for name,dx in [('L10',14),('L20',-18)]:
        shift=Affine2D().translate(dx,0)
        for artist in groups[name]:
            transform=artist.get_data_transform() if isinstance(artist,Patch) else artist.get_transform()
            artist.set_transform(transform+shift)
    return {'L10_dx_pixels':14,'L20_dx_pixels':-18}


def align_graphic_centers(figure):
    """Center the artwork below fixed headings, using shapes as anchors."""
    ax=figure.axes[0]
    figure.canvas.draw()
    renderer=figure.canvas.get_renderer()
    groups={name:[] for name in ['left','center','right']}
    anchors={name:[] for name in groups}
    arrows=[]
    headings={'Distinct\ninterfaces','L10','L20','Interaction\ngeometry'}
    for artist in [*ax.patches,*ax.lines,*ax.texts]:
        if isinstance(artist,FancyArrowPatch):
            arrows.append(artist)
            continue
        if artist in ax.texts and artist.get_text() in headings:continue
        bounds=artist.get_window_extent(renderer)
        x=(bounds.x0+bounds.x1)/2
        name='left' if x<255 else 'center' if x<690 else 'right'
        groups[name].append(artist)
        if isinstance(artist,Patch) and (name!='right' or isinstance(artist,FancyBboxPatch)):
            anchors[name].append(artist)
    assert [len(anchors[n]) for n in groups]==[30,19,2]
    assert len(arrows)==2
    before={n:Bbox.union([a.get_window_extent(renderer) for a in anchors[n]]) for n in groups}
    centers={n:(b.y0+b.y1)/2 for n,b in before.items()}
    # The already near-aligned outer groups define the content area's center.
    target=(centers['left']+centers['right'])/2
    shifts={n:target-center for n,center in centers.items()}
    for name,artists in groups.items():
        shift=Affine2D().translate(0,shifts[name])
        for artist in artists:
            transform=artist.get_data_transform() if isinstance(artist,Patch) else artist.get_transform()
            artist.set_transform(transform+shift)
    arrow_shift=target-ax.transData.transform((0,.52))[1]
    for arrow in arrows:
        arrow.set_transform(arrow.get_data_transform()+Affine2D().translate(0,arrow_shift))
    figure.canvas.draw()
    renderer=figure.canvas.get_renderer()
    for name in groups:
        after=Bbox.union([a.get_window_extent(renderer) for a in anchors[name]])
        assert np.allclose(after.get_points(),before[name].get_points()+[0,shifts[name]])
        assert np.isclose((after.y0+after.y1)/2,target)
    for arrow in arrows:
        bounds=arrow.get_window_extent(renderer)
        assert np.isclose((bounds.y0+bounds.y1)/2,target)
    all_bounds=Bbox.union([a.get_window_extent(renderer) for artists in groups.values() for a in artists]+[a.get_window_extent(renderer) for a in arrows])
    assert 0<all_bounds.y0<all_bounds.y1<original.EXPECTED_PIXELS[1]-90
    return {'center_pixels_down':original.EXPECTED_PIXELS[1]-target,
            'group_dy_pixels_down':{n:-dy for n,dy in shifts.items()},
            'arrow_dy_pixels_down':-arrow_shift}


def main():
    source=HERE/'TOC_graphic_schematic.png'
    original_hash=hashlib.sha256(source.read_bytes()).hexdigest()
    with Image.open(source) as image:
        baseline=np.array(image.convert('RGB'))
    figure=original.build_figure()
    figure.canvas.draw()
    rerender=np.asarray(figure.canvas.buffer_rgba())[:,:,:3].copy()
    assert np.array_equal(rerender,baseline),'V1 source no longer renders the selected image exactly'
    arrows=[p for p in figure.axes[0].patches if isinstance(p,FancyArrowPatch)]
    assert len(arrows)==3
    central=[p for p in arrows if p.get_zorder()==10]
    outer=[p for p in arrows if p.get_zorder()==12]
    assert len(central)==1 and len(outer)==2
    central[0].remove()
    for arrow in outer:
        # Default two-point shrink at both ends consumes almost the entire
        # short shaft at this canvas size. Keep original positions and style.
        arrow.shrinkA=0
        arrow.shrinkB=0
    shifts=center_card_contents(figure)
    figure.canvas.draw()
    previous=np.asarray(figure.canvas.buffer_rgba())[:,:,:3].copy()
    group_shifts=tighten_center_groups(figure)
    for arrow,(x0,x1) in zip(outer,[(.268,.309),(.672,.708)]):
        arrow.set_positions((x0,.52),(x1,.52))
    figure.canvas.draw()
    result=np.asarray(figure.canvas.buffer_rgba())[:,:,:3].copy()
    latest_changed=np.any(result!=previous,axis=2)
    latest_allowed=np.zeros(latest_changed.shape,dtype=bool)
    latest_allowed[45:390,255:699]=True
    assert not np.any(latest_changed & ~latest_allowed),'Spacing edit changed another panel or heading'
    before_heading=result.copy()
    headings=[t for t in figure.axes[0].texts if t.get_text()=='Local targets']
    assert len(headings)==1
    heading=headings[0]
    heading.set_text('Interaction\ngeometry')
    heading.set_position((.850,.915))
    heading.set_linespacing(.85)
    figure.canvas.draw()
    result=np.asarray(figure.canvas.buffer_rgba())[:,:,:3].copy()
    heading_changed=np.any(result!=before_heading,axis=2)
    heading_allowed=np.zeros(heading_changed.shape,dtype=bool)
    heading_allowed[0:90,700:960]=True
    assert not np.any(heading_changed & ~heading_allowed),'Heading edit changed another object'
    heading_bounds=heading.get_window_extent(figure.canvas.get_renderer())
    assert 700<heading_bounds.x0<heading_bounds.x1<960
    assert original.EXPECTED_PIXELS[1]-90<heading_bounds.y0<heading_bounds.y1<original.EXPECTED_PIXELS[1]
    changed=np.any(result!=baseline,axis=2)
    boxes=[(255,45,699,390),
           (700,107,959,276),(700,314,959,470),(700,0,960,90)]
    allowed=np.zeros(changed.shape,dtype=bool)
    for x0,y0,x1,y1 in boxes:allowed[y0:y1,x0:x1]=True
    assert not np.any(changed & ~allowed),'Unexpected change outside requested group/arrow/card regions'
    assert all(np.any(changed[y0:y1,x0:x1]) for x0,y0,x1,y1 in boxes)
    before_alignment=result.copy()
    alignment=align_graphic_centers(figure)
    result=np.asarray(figure.canvas.buffer_rgba())[:,:,:3].copy()
    assert np.array_equal(result[:90],before_alignment[:90]),'Vertical alignment changed a heading'
    before_typography=result.copy()
    panel_headings=[t for t in figure.axes[0].texts if t.get_text() in
                    {'Distinct\ninterfaces','Interaction\ngeometry'}]
    assert len(panel_headings)==2
    for title in panel_headings:
        title.set_fontsize(7.2)
        title.set_fontweight('bold')
        title.set_linespacing(.95)
        title.set_verticalalignment('center')
        title.set_position((title.get_position()[0],.915))
    figure.canvas.draw()
    result=np.asarray(figure.canvas.buffer_rgba())[:,:,:3].copy()
    typography_changed=np.any(result!=before_typography,axis=2)
    typography_allowed=np.zeros(typography_changed.shape,dtype=bool)
    typography_allowed[:90,:255]=True
    typography_allowed[:90,700:]=True
    assert not np.any(typography_changed & ~typography_allowed),'Typography edit changed another object'
    title_bounds=[t.get_window_extent(figure.canvas.get_renderer()) for t in panel_headings]
    assert np.isclose(title_bounds[0].y0,title_bounds[1].y0)
    assert np.isclose(title_bounds[0].y1,title_bounds[1].y1)
    before_label_swap=result.copy()
    middle_title=next(t for t in figure.axes[0].texts if t.get_text()=='Enzyme-specific\ninterface remodeling')
    length_labels=[t for t in figure.axes[0].texts if t.get_text() in {'L10','L20'}]
    renderer=figure.canvas.get_renderer()
    old_title_bounds=middle_title.get_window_extent(renderer)
    bottom_center=(old_title_bounds.y0+old_title_bounds.y1)/2
    top_center=(title_bounds[0].y0+title_bounds[0].y1)/2
    for text,target in [(middle_title,top_center),*[(t,bottom_center) for t in length_labels]]:
        bounds=text.get_window_extent(renderer)
        dy=target-(bounds.y0+bounds.y1)/2
        text.set_transform(text.get_transform()+Affine2D().translate(0,dy))
    figure.canvas.draw()
    result=np.asarray(figure.canvas.buffer_rgba())[:,:,:3].copy()
    swap_changed=np.any(result!=before_label_swap,axis=2)
    swap_allowed=np.zeros(swap_changed.shape,dtype=bool)
    swap_allowed[:90,295:665]=True
    swap_allowed[425:500,295:665]=True
    assert not np.any(swap_changed & ~swap_allowed),'Label swap changed artwork or another panel'
    changed=np.any(result!=baseline,axis=2)
    # Use one PNG export only. Do not call the legacy dual PNG/TIFF exporter.
    figure.savefig(OUTPUT,format='png',dpi=original.DPI,facecolor='white')
    with Image.open(OUTPUT) as final:
        assert final.size==original.EXPECTED_PIXELS
        assert np.array_equal(np.array(final.convert('RGB')),result)
        if final.mode=='RGBA':assert np.all(np.array(final)[:,:,3]==255)
    assert hashlib.sha256(source.read_bytes()).hexdigest()==original_hash
    original.plt.close(figure)
    print(json.dumps({'output':str(OUTPUT),'size':list(original.EXPECTED_PIXELS),
                      'source_sha256':original_hash,'changed_pixels':int(changed.sum()),
                      'pre_alignment_pixels_outside_requested_regions':0,'card_translations':shifts,
                      'center_group_translations':group_shifts,
                      'spacing_edit_pixels_outside_center_panel':0,
                      'right_heading':'Interaction geometry',
                      'heading_edit_pixels_outside_title_area':0,
                      'vertical_alignment':alignment,
                      'panel_heading_fontsize_pt':7.2,'panel_heading_linespacing':.95,
                      'typography_edit_pixels_outside_heading_areas':0,
                      'middle_title_above':True,'length_labels_below':True,
                      'label_swap_pixels_outside_text_areas':0,
                      'source_preserved':True,'central_arrows':0,'outer_arrows':2},indent=2))


if __name__=='__main__':
    raise SystemExit("Drawing source only. Use make_figure.py --output NEW_DIR/TOC_graphic.png")
