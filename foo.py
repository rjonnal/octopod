# Recorded script from Mayavi2
from numpy import array
try:
    engine = mayavi.engine
except NameError:
    from mayavi.api import Engine
    engine = Engine()
    engine.start()
if len(engine.scenes) == 0:
    engine.new_scene()
# ------------------------------------------- 
image_plane_widget = engine.scenes[0].children[0].children[0].children[1]
image_plane_widget.ipw.origin = array([  0.5       ,   0.5       ,  89.77268831])
image_plane_widget.ipw.slice_index = 89
image_plane_widget.ipw.slice_position = 89.77268831276514
image_plane_widget.ipw.point1 = array([ 190.5       ,    0.5       ,   89.77268831])
image_plane_widget.ipw.point2 = array([   0.5       ,  175.5       ,   89.77268831])
image_plane_widget.ipw.origin = array([  0.5       ,   0.5       ,  89.77268831])
image_plane_widget.ipw.point1 = array([ 190.5       ,    0.5       ,   89.77268831])
image_plane_widget.ipw.point2 = array([   0.5       ,  175.5       ,   89.77268831])
image_plane_widget.ipw.origin = array([  0.5       ,   0.5       ,  80.91664291])
image_plane_widget.ipw.slice_index = 80
image_plane_widget.ipw.slice_position = 80.91664290620359
image_plane_widget.ipw.point1 = array([ 190.5       ,    0.5       ,   80.91664291])
image_plane_widget.ipw.point2 = array([   0.5       ,  175.5       ,   80.91664291])
image_plane_widget.ipw.origin = array([  0.5       ,   0.5       ,  80.91664291])
image_plane_widget.ipw.point1 = array([ 190.5       ,    0.5       ,   80.91664291])
image_plane_widget.ipw.point2 = array([   0.5       ,  175.5       ,   80.91664291])
image_plane_widget.ipw.origin = array([   0.5       ,    0.5       ,  104.95448071])
image_plane_widget.ipw.slice_index = 104
image_plane_widget.ipw.slice_position = 104.95448070675768
image_plane_widget.ipw.point1 = array([ 190.5       ,    0.5       ,  104.95448071])
image_plane_widget.ipw.point2 = array([   0.5       ,  175.5       ,  104.95448071])
image_plane_widget.ipw.origin = array([   0.5       ,    0.5       ,  104.95448071])
image_plane_widget.ipw.point1 = array([ 190.5       ,    0.5       ,  104.95448071])
image_plane_widget.ipw.point2 = array([   0.5       ,  175.5       ,  104.95448071])
image_plane_widget.ipw.origin = array([   0.5,    0.5,  107. ])
image_plane_widget.ipw.slice_index = 106
image_plane_widget.ipw.slice_position = 107.0
image_plane_widget.ipw.point1 = array([ 190.5,    0.5,  107. ])
image_plane_widget.ipw.point2 = array([   0.5,  175.5,  107. ])
image_plane_widget.ipw.origin = array([   0.5,    0.5,  107. ])
image_plane_widget.ipw.point1 = array([ 190.5,    0.5,  107. ])
image_plane_widget.ipw.point2 = array([   0.5,  175.5,  107. ])
image_plane_widget.ipw.origin = array([  0.5      ,   0.5      ,  57.0190253])
image_plane_widget.ipw.slice_index = 56
image_plane_widget.ipw.slice_position = 57.019025298014846
image_plane_widget.ipw.point1 = array([ 190.5      ,    0.5      ,   57.0190253])
image_plane_widget.ipw.point2 = array([   0.5      ,  175.5      ,   57.0190253])
image_plane_widget.ipw.origin = array([  0.5      ,   0.5      ,  57.0190253])
image_plane_widget.ipw.point1 = array([ 190.5      ,    0.5      ,   57.0190253])
image_plane_widget.ipw.point2 = array([   0.5      ,  175.5      ,   57.0190253])
scene = engine.scenes[0]
scene.scene.camera.position = [474.1741250984158, 466.2043535058537, 121.24423379774402]
scene.scene.camera.focal_point = [95.5, 88.0, 54.0]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [-0.090426890690081904, -0.08587158028263514, 0.99219405820630002]
scene.scene.camera.clipping_range = [266.05290475593517, 884.79304367314751]
scene.scene.camera.compute_view_plane_normal()
scene.scene.render()
scene.scene.camera.position = [331.03321626867324, 352.43551533719966, 460.88127908686238]
scene.scene.camera.focal_point = [95.5, 88.0, 54.0]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [-0.51944320709572334, -0.54719589076029085, 0.656319595728454]
scene.scene.camera.clipping_range = [286.53720140256542, 859.00753457776682]
scene.scene.camera.compute_view_plane_normal()
scene.scene.render()
