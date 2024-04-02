# pylint: disable=no-name-in-module
"""
Defines 3D models.
"""

import os
import tempfile
import warnings
import webbrowser
from datetime import datetime
from typing import List

try:
    import gmsh
except (TypeError, OSError):
    pass
import matplotlib.pyplot as plt
import numpy as np
# pylint: disable=no-name-in-module
from OCP.Interface import Interface_Static
import OCP.IFSelect
from OCP.STEPControl import STEPControl_Reader
from OCP.TopoDS import TopoDS_Shape

import dessia_common.core as dc
from dessia_common.errors import ConsistencyError
import dessia_common.files as dcf
import volmdlr
import volmdlr.templates
from volmdlr.utils.step_writer import STEP_HEADER, STEP_FOOTER
from volmdlr import from_ocp
from volmdlr.core import Primitive3D, BoundingBox, get_babylon_data
from volmdlr.shells import Shell3D


class VolumeModel(dc.PhysicalObject):
    """
    A class containing one or several :class:`volmdlr.core.Primitive3D`.

    :param primitives: The vector's abscissa
    :type primitives: List[:class:`volmdlr.core.Primitive3D`]
    :param name: The VolumeModel's name
    :type name: str
    """
    _standalone_in_db = True
    _eq_is_data_eq = True
    _non_serializable_attributes = ['shells', 'bounding_box']
    _non_data_eq_attributes = ['name', 'shells', 'bounding_box', 'contours',
                               'faces']
    _non_data_hash_attributes = ['name', 'shells', 'bounding_box', 'contours',
                                 'faces']
    _dessia_methods = ['to_stl_model']

    def __init__(self, primitives: List[Primitive3D], name: str = ''):
        self.primitives = primitives
        self.name = name
        self.shells = []
        self._bbox = None
        dc.PhysicalObject.__init__(self, name=name)

    def __hash__(self):
        return sum(hash(point) for point in self.primitives)

    def __eq__(self, other):
        if self.__class__.__name__ != other.__class__.__name__:
            return False
        equ = True
        if len(self.primitives) != len(other.primitives):
            return False
        for p1, p2 in zip(self.primitives, other.primitives):
            # TODO: if 2 volume models has exact same primitives but in a different order, they are different?
            equ = equ and p1 == p2
        return equ

    @property
    def bounding_box(self):
        """
        Returns the bounding box.

        """
        if not self._bbox:
            self._bbox = self._bounding_box()
        return self._bbox

    @bounding_box.setter
    def bounding_box(self, new_bounding_box):
        self._bbox = new_bounding_box

    def _bounding_box(self) -> BoundingBox:
        """
        Computes the bounding box of the model.
        """
        return BoundingBox.from_bounding_boxes([p.bounding_box for p in self.primitives if hasattr(p, "bounding_box")])

    def volume(self) -> float:
        """
        Return the sum of volumes of the primitives.

        It does not make any Boolean operation in case of overlapping.

        """
        volume = 0
        for primitive in self.primitives:
            volume += primitive.volume()
        return volume

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D,
                 angle: float):
        """
        Rotates the VolumeModel.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: angle rotation
        :return: a new rotated VolumeModel
        """
        new_primitives = [
            primitive.rotation(center, axis, angle) for
            primitive in self.primitives]
        return VolumeModel(new_primitives, self.name)

    def translation(self, offset: volmdlr.Vector3D):
        """
        Translates the VolumeModel.

        :param offset: translation vector
        :return: A new translated VolumeModel
        """
        new_primitives = [primitive.translation(offset) for
                          primitive in self.primitives]
        return VolumeModel(new_primitives, self.name)

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new VolumeModel.

        side = 'old' or 'new'
        """
        new_primitives = [primitive.frame_mapping(frame, side)
                          for primitive in self.primitives]
        return VolumeModel(new_primitives, self.name)

    def copy(self, deep=True, memo=None):
        """
        Specific copy.
        """
        new_primitives = [primitive.copy(deep=deep, memo=memo) for primitive in self.primitives]
        return VolumeModel(new_primitives, self.name)

    def plot(self, equal_aspect=True):
        """
        Matplotlib plot of model.

        To use for debug.
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d', adjustable='box')
        for primitive in self.primitives:
            primitive.plot(ax)
        if not equal_aspect:
            # ax.set_aspect('equal')
            ax.set_aspect('auto')
        ax.margins(0.1)
        return ax

    def babylon_data(self, merge_meshes=True):
        """
        Get babylonjs data.

        :return: Dictionary with babylon data.
        """

        return get_babylon_data(self, merge_meshes=merge_meshes)

    @classmethod
    def babylonjs_script(cls, babylon_data, use_cdn=True, **kwargs):
        """
        Run babylonjs script.

        """
        if use_cdn:
            script = volmdlr.templates.BABYLON_UNPACKER_CDN_HEADER  # .substitute(name=page_name)
        else:
            script = volmdlr.templates.BABYLON_UNPACKER_EMBEDDED_HEADER  # .substitute(name=page_name)

        script += volmdlr.templates.BABYLON_UNPACKER_BODY_TEMPLATE.substitute(
            babylon_data=babylon_data)
        return script

    def babylonjs(
            self,
            page_name: str = None,
            use_cdn: bool = True,
            debug: bool = False,
            merge_meshes: bool = True,
            dark_mode: bool = False,
    ):
        """
        Generate and display an HTML file to visualize the 3D model using Babylon.js in a web browser.

        This method creates a 3D representation of the volume model using the Babylon.js framework.
        The method allows options for debugging, merging meshes, and toggling dark mode for the visualization.
        The resulting HTML file can either be a temporary file or a user-specified file.

        :param page_name: The name of the HTML file to be generated. If None, a temporary file is created.
        :type page_name: str, optional
        :param use_cdn: Flag to use CDN for loading Babylon.js resources. Defaults to True.
        :type use_cdn: bool
        :param debug: Enable debugging mode for more detailed console output in the browser. Defaults to False.
        :type debug: bool
        :param merge_meshes: Flag to chose to merge all the faces of each shell into a single mesh. Defaults to True.
            If False, shell are decomposed according to their faces in the Babylon.js scene nodes tree.
        :type merge_meshes: bool
        :param dark_mode: Enable dark mode for the HTML visualization. Defaults to False.
        :type dark_mode: bool

        :return: The file path of the generated HTML file.
        :rtype: str
        """
        babylon_data = self.babylon_data(merge_meshes=merge_meshes)
        babylon_data["dark_mode"] = 1 if dark_mode else 0
        script = self.babylonjs_script(babylon_data, use_cdn=use_cdn, debug=debug)
        if page_name is None:
            with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as file:
                file.write(bytes(script, "utf8"))
            page_name = file.name
        else:
            if not page_name.endswith(".html"):
                page_name += ".html"
            with open(page_name, "w", encoding="utf-8") as file:
                file.write(script)

        webbrowser.open("file://" + os.path.realpath(page_name))

        return page_name

    def save_babylonjs_to_file(self, filename: str = None, use_cdn=True, debug=False, dark_mode=False):
        """Export a html file of the model."""
        babylon_data = self.babylon_data()
        babylon_data['dark_mode'] = 1 if dark_mode else 0
        script = self.babylonjs_script(babylon_data, use_cdn=use_cdn, debug=debug)
        if filename is None:
            with tempfile.NamedTemporaryFile(suffix=".html",
                                             delete=False) as file:
                file.write(bytes(script, 'utf8'))
                return file.name

        if not filename.endswith('.html'):
            filename += '.html'

        with open(filename, 'w', encoding='utf-8') as file:
            file.write(script)
        return filename

    def to_mesh_list(self):
        """
        Converts the volume model to a list Mesh3D object.

        :return: A list of Mesh3D objects representing the VolumeModel shells.
        """
        meshes = []

        for shell in self.get_shells():
            mesh = shell.triangulation()

            if len(mesh.triangles) > 0:
                meshes.append(mesh)
                meshes[-1].name = shell.name

        return meshes

    def to_mesh(self, merge_vertices: bool = True, merge_triangles: bool = True):
        """
        Converts the volume model to a Mesh3D object.

        :param merge_vertices: Flag to indicate whether to merge vertices of the shells meshes.
        :param merge_triangles: Flag to indicate whether to merge triangles of the shells meshes.

        :return: A Mesh3D of the VolumeModel
        """
        meshes = self.to_mesh_list()

        if not meshes:
            raise ValueError("VolumeModel has no primitive that can be converted to mesh.")

        merged_mesh = meshes[0]
        for mesh in meshes[1:]:
            merged_mesh = merged_mesh.merge(mesh, merge_vertices=merge_vertices, merge_triangles=merge_triangles)

        merged_mesh.name = self.name

        return merged_mesh

    def to_stl_model(self):
        """Converts the model into a stl object."""
        warnings.warn(
            "volmdlr.stl module is deprecated. Use volmdlr.display module and 'Mesh3D' class instead for STL export.",
            DeprecationWarning
        )

        mesh = self.to_mesh()

        # from volmdlr import stl
        stl = volmdlr.stl.Stl.from_display_mesh(mesh)
        return stl

    def to_stl(self, filepath: str):
        """Export a stl file of the model."""
        self.to_mesh().save_to_stl_file(filepath)

    def to_stl_stream(self, stream: dcf.BinaryFile):
        """Converts the model into a stl stream file."""
        self.to_mesh().save_to_stl_stream(stream)
        return stream

    def to_step(self, filepath: str):
        """Export a step file of the model."""
        if not (filepath.endswith('.step') or filepath.endswith('.stp')):
            filepath += '.step'
        with open(filepath, 'w', encoding='utf-8') as file:
            self.to_step_stream(file)

    def to_step_stream(self, stream: dcf.StringFile):
        """
        Export object CAD to given stream in STEP format.

        """
        step_content = STEP_HEADER.format(name=self.name,
                                          filename='',
                                          timestamp=datetime.now().isoformat(),
                                          version=volmdlr.__version__)
        current_id = 2

        for primitive in self.primitives:
            if primitive.__class__.__name__ in ('OpenShell3D', 'ClosedShell3D') or hasattr(primitive, "shell_faces"):
                primitive_content, primitive_id, _ = primitive.to_step_product(current_id)
            elif primitive.__class__.__name__ in ('Assembly', 'Compound'):
                primitive_content, primitive_id, _ = primitive.to_step(current_id)
            else:
                continue

            step_content += primitive_content
            current_id = primitive_id

        step_content += STEP_FOOTER

        stream.write(step_content)

    @classmethod
    def from_step(cls, step_file: str, name: str = ""):
        """
        Translates a STEP file into a volume model using Open CASCADE Technology (OCCT).

        This method reads a STEP file, converts it into a series of shapes using OCCT, and then translates these shapes
         into a volume model. The unit of measurement is set to meters.

        Note: This method is experimental and may contain bugs.
        """
        # Set the unit to meter
        reader = STEPControl_Reader()
        Interface_Static.SetCVal_s("xstep.cascade.unit", "M")
        read_status = reader.ReadFile(step_file)
        if read_status != OCP.IFSelect.IFSelect_RetDone:
            raise ValueError("STEP File could not be loaded")
        for i in range(reader.NbRootsForTransfer()):
            reader.TransferRoot(i + 1)

        occ_shapes = []
        for i in range(reader.NbShapes()):
            occ_shapes.append(reader.Shape(i + 1))
        return cls.from_ocp(occ_shapes, name=name)


    @classmethod
    def from_ocp(cls, ocp_shapes:List[TopoDS_Shape], name: str = ""):
        """
        Instanciate a volume model from a list of OCCT shapes.
        """
        # Make sure that we extract all the solids
        solids = []
        for shape in ocp_shapes:
            shape_type = from_ocp.shapetype(shape)
            # TODO: we get only the shells inside the Compound because circular imports
            if shape_type in (0, 1, 2):
                list_of_shells = from_ocp.get_shells(shape)
                for shell in list_of_shells:
                    solids.append(Shell3D.from_ocp(shell))
            elif from_ocp.shapetype(shape) == 3:
                Shell3D.from_ocp(shape)

        return cls(primitives=solids, name=name)

    def volmdlr_volume_model(self):
        """
        Method needed due to PhysicalObject inheritance.
        """
        return self

    def get_geo_lines(self):
        """
        Gets the lines that define a VolumeModel geometry in a .geo file.

        :return: A list of lines that describe the geometry
        :rtype: List[str]

        """

        update_data = {'point_account': 0,
                       'line_account': 0,
                       'line_loop_account': 0,
                       'surface_account': 0,
                       'surface_loop_account': 0}

        lines = []
        volume = 0
        for primitive in self.primitives:
            if isinstance(primitive, volmdlr.shells.ClosedShell3D):
                volume += 1
                lines_primitives, update_data = primitive.get_geo_lines(update_data)
                lines.extend(lines_primitives)
                surface_loop = ((lines[-1].split('('))[1].split(')')[0])
                lines.append('Volume(' + str(volume) + ') = {' + surface_loop + '};')
            elif isinstance(primitive, volmdlr.shells.OpenShell3D):
                lines_primitives, update_data = primitive.get_geo_lines(update_data)
                lines.extend(lines_primitives)

        return lines

    def get_mesh_lines(self, factor: float, **kwargs):
        """
        Gets the lines that define mesh parameters for a VolumeModel, to be added to a .geo file.

        :param factor: A float, between 0 and 1, that describes the mesh quality
        (1 for coarse mesh - 0 for fine mesh)
        :type factor: float
        :param curvature_mesh_size: Activate the calculation of mesh element sizes based on curvature
        (with curvature_mesh_size elements per 2*Pi radians), defaults to 0
        :type curvature_mesh_size: int, optional
        :param min_points: Check if there are enough points on small edges (if it is not, we force to have min_points
        on that edge), defaults to None
        :type min_points: int, optional
        :param initial_mesh_size: If factor=1, it will be initial_mesh_size elements per dimension, defaults to 5
        :type initial_mesh_size: float, optional

        :return: A list of lines that describe mesh parameters
        :rtype: List[str]
        """

        for element in [('curvature_mesh_size', 0), ('min_points', None), ('initial_mesh_size', 5)]:
            if element[0] not in kwargs:
                kwargs[element[0]] = element[1]

        field_num = 1
        field_nums = []
        lines = []

        lines.append('Mesh.CharacteristicLengthMin = 0;')
        lines.append('Mesh.CharacteristicLengthMax = 1e+22;')
        lines.append('Geometry.Tolerance = 1e-5;')
        lines.append('Mesh.AngleToleranceFacetOverlap = 0.01;')
        lines.append('General.Verbosity = 0;')

        for i, primitive in enumerate(self.primitives):
            if isinstance(primitive, volmdlr.shells.ClosedShell3D):
                bbx = primitive.bounding_box
                dim1, dim2, dim3 = (bbx.xmax - bbx.xmin), (bbx.ymax - bbx.ymin), (bbx.zmax - bbx.zmin)
                volume = dim1 * dim2 * dim3

                if factor == 0:
                    factor = 1e-3

                size = ((volume ** (1. / 3.)) / kwargs['initial_mesh_size']) * factor

                if kwargs['min_points']:
                    lines.extend(primitive.get_mesh_lines_with_transfinite_curves(min_points=kwargs['min_points'],
                                                                                  size=size))

                lines.append('Field[' + str(field_num) + '] = MathEval;')
                lines.append('Field[' + str(field_num) + '].F = "' + str(size) + '";')

                lines.append('Field[' + str(field_num + 1) + '] = Restrict;')
                lines.append('Field[' + str(field_num + 1) + '].InField = ' + str(field_num) + ';')
                lines.append('Field[' + str(field_num + 1) + '].VolumesList = {' + str(i + 1) + '};')
                field_nums.append(field_num + 1)
                field_num += 2

            elif isinstance(primitive, volmdlr.shells.OpenShell3D):
                continue

        lines.append('Field[' + str(field_num) + '] = MinAniso;')
        lines.append('Field[' + str(field_num) + '].FieldsList = {' + str(field_nums)[1:-1] + '};')
        lines.append('Background Field = ' + str(field_num) + ';')

        lines.append('Mesh.MeshSizeFromCurvature = ' + str(kwargs['curvature_mesh_size']) + ';')

        lines.append('Coherence;')

        return lines

    def to_geo_stream(self, stream: dcf.StringFile,
                      factor: float, **kwargs):
        """
        Gets the .geo file for the VolumeModel.

        :param file_name: The geo. file name
        :type file_name: str
        :param factor: A float, between 0 and 1, that describes the mesh quality
        (1 for coarse mesh - 0 for fine mesh)
        :type factor: float
        :param curvature_mesh_size: Activate the calculation of mesh element sizes based on curvature
        (with curvature_mesh_size elements per 2*Pi radians), defaults to 0
        :type curvature_mesh_size: int, optional
        :param min_points: Check if there are enough points on small edges (if it is not, we force to have min_points
        on that edge), defaults to None
        :type min_points: int, optional
        :param initial_mesh_size: If factor=1, it will be initial_mesh_size elements per dimension, defaults to 5
        :type initial_mesh_size: float, optional

        :return: A txt file
        :rtype: .txt
        """

        for element in [('curvature_mesh_size', 0), ('min_points', None), ('initial_mesh_size', 5)]:
            if element[0] not in kwargs:
                kwargs[element[0]] = element[1]

        lines = self.get_geo_lines()
        lines.extend(self.get_mesh_lines(factor,
                                         curvature_mesh_size=kwargs['curvature_mesh_size'],
                                         min_points=kwargs['min_points'],
                                         initial_mesh_size=kwargs['initial_mesh_size']))

        content = ''
        for line in lines:
            content += line + '\n'

        stream.write(content)

    def to_geo(self, file_name: str = '',
               factor: float = 0.5, **kwargs):
        # curvature_mesh_size: int = 0,
        # min_points: int = None,
        # initial_mesh_size: float = 5):
        """
        Gets the .geo file for the VolumeModel.

        :param file_name: The geo. file name
        :type file_name: str
        :param factor: A float, between 0 and 1, that describes the mesh quality
        (1 for coarse mesh - 0 for fine mesh)
        :type factor: float
        :param curvature_mesh_size: Activate the calculation of mesh element sizes based on curvature
        (with curvature_mesh_size elements per 2*Pi radians), defaults to 0
        :type curvature_mesh_size: int, optional
        :param min_points: Check if there are enough points on small edges (if it is not, we force to have min_points
        on that edge), defaults to None
        :type min_points: int, optional
        :param initial_mesh_size: If factor=1, it will be initial_mesh_size elements per dimension, defaults to 5
        :type initial_mesh_size: float, optional

        :return: A txt file
        :rtype: .txt
        """

        for element in [('curvature_mesh_size', 0), ('min_points', None), ('initial_mesh_size', 5)]:
            if element[0] not in kwargs:
                kwargs[element[0]] = element[1]

        if not (file_name.endswith('.geo') or file_name.endswith('.geo')):
            file_name += '.geo'

        with open(file_name, mode='w', encoding='utf-8') as file:

            self.to_geo_stream(file, factor,
                               curvature_mesh_size=kwargs['curvature_mesh_size'],
                               min_points=kwargs['min_points'],
                               initial_mesh_size=kwargs['initial_mesh_size'])

        # for element in [('curvature_mesh_size', 0), ('min_points', None), ('initial_mesh_size', 5)]:
        #     if element[0] not in kwargs:
        #         kwargs[element[0]] = element[1]

        # # try:
        # #     curvature_mesh_size = kwargs['curvature_mesh_size']
        # # except KeyError:
        # #     curvature_mesh_size = 0
        # # try:
        # #     min_points = kwargs['min_points']
        # # except KeyError:
        # #     min_points = None
        # # try:
        # #     initial_mesh_size = kwargs['initial_mesh_size']
        # # except KeyError:
        # #     initial_mesh_size = 5

        # lines = self.get_geo_lines()
        # lines.extend(self.get_mesh_lines(factor,
        #                                   curvature_mesh_size=kwargs['curvature_mesh_size'],
        #                                   min_points=kwargs['min_points'],
        #                                   initial_mesh_size=kwargs['initial_mesh_size']))
        # with open(file_name + '.geo', 'w', encoding="utf-8") as file:
        #     for line in lines:
        #         file.write(line)
        #         file.write('\n')
        # file.close()

    def to_geo_with_stl(self, file_name: str,
                        factor: float, **kwargs):
        # curvature_mesh_size: int = 0,
        # min_points: int = None,
        # initial_mesh_size: float = 5):
        """
        Gets the .geo file for the VolumeModel, with saving each closed shell in a stl file.

        :param file_name: The geo. file name
        :type file_name: str
        :param factor: A float, between 0 and 1, that describes the mesh quality
        (1 for coarse mesh - 0 for fine mesh)
        :type factor: float
        :param curvature_mesh_size: Activate the calculation of mesh element sizes based on curvature
        (with curvature_mesh_size elements per 2*Pi radians), defaults to 0
        :type curvature_mesh_size: int, optional
        :param min_points: Check if there are enough points on small edges (if it is not, we force to have min_points
        on that edge), defaults to None
        :type min_points: int, optional
        :param initial_mesh_size: If factor=1, it will be initial_mesh_size elements per dimension, defaults to 5
        :type initial_mesh_size: float, optional

        :return: A txt file
        :rtype: .txt
        """

        for element in [('curvature_mesh_size', 0), ('min_points', None), ('initial_mesh_size', 5)]:
            if element[0] not in kwargs:
                kwargs[element[0]] = element[1]

        # try:
        #     curvature_mesh_size = kwargs['curvature_mesh_size']
        # except KeyError:
        #     curvature_mesh_size = 0
        # try:
        #     min_points = kwargs['min_points']
        # except KeyError:
        #     min_points = None
        # try:
        #     initial_mesh_size = kwargs['initial_mesh_size']
        # except KeyError:
        #     initial_mesh_size = 5

        lines = self.get_geo_lines()
        lines.extend(self.get_mesh_lines(factor,
                                         curvature_mesh_size=kwargs['curvature_mesh_size'],
                                         min_points=kwargs['min_points'],
                                         initial_mesh_size=kwargs['initial_mesh_size']))

        contours, faces_account = [], 0
        surfaces = []
        for i, primitive in enumerate(self.primitives):
            if i == 0:
                surfaces.append(list(range(1, 1 + len(primitive.faces))))
                face_contours = [face.outer_contour3d for face in primitive.faces]
                contours.append(face_contours)
                lines.append('Mesh 2;')
                lines.append('Physical Surface(' + str(i + 1) + ') = {' + str(surfaces[i])[1:-1] + '};')
                lines.append('Save "' + file_name + '.stl" ;')
                faces_account += len(primitive.faces) + 1
            else:
                surfaces.append(list(range(faces_account, faces_account + len(primitive.faces))))
                face_contours = [face.outer_contour3d for face in primitive.faces]
                surfaces = self.update_surfaces_list(face_contours, surfaces, contours, i)
                # for k, face_c in enumerate(face_contours):
                #     for l, contour_l in enumerate(contours):
                #         for c, contour in enumerate(contour_l):
                #             if face_c.is_superposing(contour):
                #                 surfaces[i][k] = surfaces[l][c]
                #                 continue
                lines.append('Mesh 2;')
                lines.append('Physical Surface(' + str(i + 1) + ') = {' + str(surfaces[i])[1:-1] + '};')
                lines.append('Save "' + file_name + '.stl" ;')
                faces_account += len(primitive.faces) + 1
                contours.append(face_contours)

        return lines

    @staticmethod
    def update_surfaces_list(face_contours, surfaces, contours, i):
        """Update surfaces list."""
        for k_f, face_c in enumerate(face_contours):
            for l_c, contour_l in enumerate(contours):
                for c_c, contour in enumerate(contour_l):
                    if face_c.is_superposing(contour):
                        surfaces[i][k_f] = surfaces[l_c][c_c]
                        continue
        return surfaces

    def to_msh(self, mesh_dimension: int, factor: float,
               mesh_order: int = 1, file_name: str = '', **kwargs):
        # curvature_mesh_size: int = 0,
        # min_points: int = None,
        # initial_mesh_size: float = 5):
        """
        Gets .msh file for the VolumeModel generated by gmsh.

        :param file_name: The msh. file name
        :type file_name: str
        :param mesh_dimension: The mesh dimension (1: 1D-Edge, 2: 2D-Triangle, 3D-Tetrahedra)
        :type mesh_dimension: int
        :param mesh_order: The msh order (1: linear, 2: 2nd order)
        :type mesh_order: int
        :param factor: A float, between 0 and 1, that describes the mesh quality
        (1 for coarse mesh - 0 for fine mesh)
        :type factor: float
        :param curvature_mesh_size: Activate the calculation of mesh element sizes based on curvature
        (with curvature_mesh_size elements per 2*Pi radians), defaults to 0
        :type curvature_mesh_size: int, optional
        :param min_points: Check if there are enough points on small edges (if it is not, we force to have min_points
        on that edge), defaults to None
        :type min_points: int, optional
        :param initial_mesh_size: If factor=1, it will be initial_mesh_size elements per dimension, defaults to 5
        :type initial_mesh_size: float, optional

        :return: A txt file
        :rtype: .txt
        """

        for element in [('curvature_mesh_size', 0), ('min_points', None), ('initial_mesh_size', 5)]:
            if element[0] not in kwargs:
                kwargs[element[0]] = element[1]

        # try:
        #     curvature_mesh_size = kwargs['curvature_mesh_size']
        # except KeyError:
        #     curvature_mesh_size = 0
        # try:
        #     min_points = kwargs['min_points']
        # except KeyError:
        #     min_points = None
        # try:
        #     initial_mesh_size = kwargs['initial_mesh_size']
        # except KeyError:
        #     initial_mesh_size = 5

        if file_name == '':
            with tempfile.NamedTemporaryFile(delete=False) as file:
                file_name = file.name

        self.to_geo(file_name=file_name,
                    factor=factor,
                    curvature_mesh_size=kwargs['curvature_mesh_size'],
                    min_points=kwargs['min_points'],
                    initial_mesh_size=kwargs['initial_mesh_size'])

        self.generate_msh_file(file_name, mesh_dimension, mesh_order)

        # gmsh.initialize()
        # gmsh.open(file_name + ".geo")

        # gmsh.model.geo.synchronize()
        # gmsh.model.mesh.generate(mesh_dimension)

        # gmsh.write(file_name + ".msh")

        # gmsh.finalize()

    @staticmethod
    def generate_msh_file(file_name, mesh_dimension, mesh_order):
        """
        Generates a mesh written in a .msh file using GMSH library.

        :param file_name: DESCRIPTION
        :type file_name: TYPE
        :param mesh_dimension: DESCRIPTION
        :type mesh_dimension: TYPE
        :param mesh_order: The msh order (1: linear, 2: 2nd order)
        :type mesh_order: int

        :return: DESCRIPTION
        :rtype: TYPE

        """

        gmsh.initialize()
        gmsh.open(file_name + ".geo")

        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(mesh_dimension)
        gmsh.model.mesh.setOrder(mesh_order)

        gmsh.write(file_name + ".msh")

        gmsh.finalize()

    def to_msh_stream(self, mesh_dimension: int,
                      factor: float, stream: dcf.StringFile,
                      mesh_order: int = 1,
                      file_name: str = '', **kwargs):
        """
        Gets .msh file for the VolumeModel generated by gmsh.

        :param file_name: The msh. file name
        :type file_name: str
        :param mesh_dimension: The mesh dimension (1: 1D-Edge, 2: 2D-Triangle, 3D-Tetrahedra)
        :type mesh_dimension: int
        :param mesh_order: The msh order (1: linear, 2: 2nd order)
        :type mesh_order: int

        :param factor: A float, between 0 and 1, that describes the mesh quality
        (1 for coarse mesh - 0 for fine mesh)
        :type factor: float
        :param curvature_mesh_size: Activate the calculation of mesh element sizes based on curvature
        (with curvature_mesh_size elements per 2*Pi radians), defaults to 0
        :type curvature_mesh_size: int, optional
        :param min_points: Check if there are enough points on small edges (if it is not, we force to have min_points
        on that edge), defaults to None
        :type min_points: int, optional
        :param initial_mesh_size: If factor=1, it will be initial_mesh_size elements per dimension, defaults to 5
        :type initial_mesh_size: float, optional

        :return: A txt file
        :rtype: .txt
        """

        for element in [('curvature_mesh_size', 0), ('min_points', None), ('initial_mesh_size', 5)]:
            if element[0] not in kwargs:
                kwargs[element[0]] = element[1]

        if file_name == '':
            with tempfile.NamedTemporaryFile(delete=False) as file:
                file_name = file.name

        self.to_geo(file_name=file_name,
                    factor=factor,
                    curvature_mesh_size=kwargs['curvature_mesh_size'],
                    min_points=kwargs['min_points'],
                    initial_mesh_size=kwargs['initial_mesh_size'])

        gmsh.initialize()
        gmsh.open(file_name + ".geo")

        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(mesh_dimension)
        gmsh.model.mesh.setOrder(mesh_order)

        lines = []
        lines.append('$MeshFormat')
        lines.append('4.1 0 8')
        lines.append('$EndMeshFormat')

        lines.extend(self.get_nodes_lines(gmsh))
        lines.extend(self.get_elements_lines(gmsh))

        content = ''
        for line in lines:
            content += line + '\n'

        stream.write(content)

        # gmsh.finalize()

    def to_msh_file(self, mesh_dimension: int,
                    factor: float, stream: dcf.StringFile,
                    mesh_order: int = 1, file_name: str = '', **kwargs):
        """ Convert and write model to a .msh file. """

        for element in [('curvature_mesh_size', 0), ('min_points', None), ('initial_mesh_size', 5)]:
            if element[0] not in kwargs:
                kwargs[element[0]] = element[1]

        if file_name == '':
            with tempfile.NamedTemporaryFile(delete=False) as file:
                file_name = file.name

        with open(file_name, mode='w', encoding='utf-8') as file:
            self.to_msh_stream(mesh_dimension=mesh_dimension,
                               factor=factor, file_name=file_name,
                               mesh_order=mesh_order,
                               stream=stream,
                               curvature_mesh_size=kwargs['curvature_mesh_size'],
                               min_points=kwargs['min_points'],
                               initial_mesh_size=kwargs['initial_mesh_size'])

    @staticmethod
    def get_nodes_lines(gmsh_model):
        """Get nodes lines."""
        lines_nodes = []
        lines_nodes.append('$Nodes')

        tag = None
        entities = gmsh_model.model.getEntities()
        for dim, tag in entities:
            node_tags, node_coords, _ = gmsh_model.model.mesh.getNodes(dim, tag)

            lines_nodes.append(str(dim) + ' ' + str(tag) + ' ' + '0 ' + str(len(node_tags)))
            for tag in node_tags:
                lines_nodes.append(str(tag))
            for n in range(0, len(node_coords), 3):
                lines_nodes.append(str(node_coords[n:n + 3])[1:-1])

        lines_nodes.insert(1, str(len(entities)) + ' ' + str(tag) + ' 1 ' + str(tag))
        lines_nodes.append('$EndNodes')

        return lines_nodes

    @staticmethod
    def get_elements_lines(gmsh_model):
        """
        Helper function to export the volume model into gmsh format.
        """
        lines_elements = []
        lines_elements.append('$Elements')

        entities = gmsh_model.model.getEntities()
        if entities:
            for dim, tag in entities:
                elem_types, elem_tags, elem_node_tags = gmsh_model.model.mesh.getElements(dim, tag)

                lines_elements.append(str(dim) + ' ' + str(tag) + ' ' + str(elem_types[0]) +
                                      ' ' + str(len(elem_tags[0])))
                range_list = int(len(elem_node_tags[0]) / len(elem_tags[0]))
                for n in range(0, len(elem_node_tags[0]), range_list):
                    lines_elements.append(str(elem_tags[0][int(n / range_list)]) + ' ' +
                                          str(elem_node_tags[0][n:n + range_list])[1:-1])

            tag = str(elem_tags[0][int(n / range_list)])
            lines_elements.insert(1, str(len(entities)) + ' ' + tag + ' 1 ' + tag)
        lines_elements.append('$EndElements')

        return lines_elements

    def get_shells(self):
        """
        Dissociates all the assemblies to get a list of shells only.

        :return: A list of closed shells
        :rtype: List[OpenShell3D]
        """

        list_shells = []

        def unpack_assembly(assembly):
            for prim in assembly.primitives:
                if prim.__class__.__name__ in ('Assembly', "Compound"):
                    unpack_assembly(prim)
                elif hasattr(prim, "faces") or hasattr(prim, "shell_faces"):
                    list_shells.append(prim)

        for primitive in self.primitives:
            if primitive.__class__.__name__ in ('Assembly', "Compound"):
                unpack_assembly(primitive)
            elif hasattr(primitive, "faces") or hasattr(primitive, "shell_faces"):
                list_shells.append(primitive)

        return list_shells


class MovingVolumeModel(VolumeModel):
    """
    A volume model with possibility to declare time steps at which the primitives are positioned with frames.

    """

    def __init__(self, primitives: List[Primitive3D], step_frames: List[List[volmdlr.Frame3D]], name: str = ''):
        VolumeModel.__init__(self, primitives=primitives, name=name)
        self.step_frames = step_frames

        if not self.is_consistent():
            raise ConsistencyError

    def is_consistent(self):
        """ Check if the number of frames for each step corresponds to the number of primitives of the model. """
        n_primitives = len(self.primitives)
        for frames in self.step_frames:
            if len(frames) != n_primitives:
                return False
        return True

    def step_volume_model(self, istep: int):
        """
        Moves the volume model with a list of local frames.
        """
        primitives = []
        for primitive, frame in zip(self.primitives, self.step_frames[istep]):
            primitives.append(
                primitive.frame_mapping(frame, side='old'))
        return VolumeModel(primitives)

    def babylon_data(self, merge_meshes=True):
        """
        Get babylonjs data.

        :return: Dictionary with babylonjs data.
        """
        meshes = []
        primitives_to_meshes = []
        for i_prim, primitive in enumerate(self.primitives):
            if hasattr(primitive, 'babylon_meshes'):
                meshes.extend(primitive.babylon_meshes(merge_meshes=merge_meshes))
                primitives_to_meshes.append(i_prim)

        # Compute max length in each direction
        all_positions = []
        for mesh in meshes:
            positions = mesh["positions"]
            all_positions.extend(positions)

        # Convert to a NumPy array and reshape
        positions_array = np.array(all_positions).reshape(-1, 3)

        # Compute min and max for each dimension
        min_vals = positions_array.min(axis=0)
        max_vals = positions_array.max(axis=0)

        # Calculate max length of the bounding box
        max_length = np.max(max_vals - min_vals)

        # Calculate center point of the bounding box
        center = (0.5 * (min_vals + max_vals)).tolist()

        steps = []
        for istep, frames in enumerate(self.step_frames):

            # step_positions = []
            # step_orientations = []
            step = {'time': istep}
            for iframe, frame in enumerate(frames):
                if iframe in primitives_to_meshes:
                    imesh = primitives_to_meshes.index(iframe)
                    step[imesh] = {}
                    step[imesh]['position'] = list(round(frame.origin, 6))
                    step[imesh]['orientations'] = [list(round(frame.u, 6)),
                                                   list(round(frame.v, 6)),
                                                   list(round(frame.w, 6))]

            steps.append(step)

        babylon_data = {'meshes': meshes,
                        'max_length': max_length,
                        'center': center,
                        'steps': steps}
        return babylon_data
