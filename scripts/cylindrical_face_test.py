import volmdlr
import volmdlr.faces
import volmdlr.wires as vme
import volmdlr.step

surface = volmdlr.faces.ToroidalSurface3D.load_from_file(r'C:\Users\gabri\Documents\dessia\GitHub\volmdlr\scripts\toroidal_surf_27.json')
contour = vme.Contour3D.load_from_file(r'C:\Users\gabri\Documents\dessia\GitHub\volmdlr\scripts\contour_surf_27.json')

face = surface.face_from_contours3d([contour])
ax = face.surface2d.plot()
ax.set_aspect('auto')
face.outer_contour3d.plot()
face.babylonjs()

