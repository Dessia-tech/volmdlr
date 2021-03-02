
import os
import pkg_resources
from string import Template


babylon_unpacker_cdn_header = '''
<!doctype html>
<html>
<head>
   <meta charset="utf-8">
   <title>Babylon from volmdlr</title>
   <style>
      html, body {
         overflow: hidden;
         width: 100%;
         height: 100%;
         margin: 0;
         padding: 0;
      }
      #renderCanvas {
         width: 100%;
         height: 100%;
         touch-action: none;
      }
   </style>
      <script src="https://preview.babylonjs.com/babylon.js"></script>
      <script src="https://preview.babylonjs.com/loaders/babylonjs.loaders.min.js"></script>
      <script src="https://code.jquery.com/pep/0.4.3/pep.js"></script>
      <script src='https://unpkg.com/earcut@2.1.1/dist/earcut.min.js'></script>
      <script src='https://preview.babylonjs.com/gui/babylon.gui.min.js'></script>
</head>
'''

babylon_unpacker_embedded_header = '''
<!doctype html>
<html>
<head>
   <meta charset="utf-8">
   <title>Babylon from volmdlr</title>
   <style>
      html, body {
         overflow: hidden;
         width: 100%;
         height: 100%;
         margin: 0;
         padding: 0;
      }
      #renderCanvas {
         width: 100%;
         height: 100%;
         touch-action: none;
      }
   </style>
   <script>
   '''

for filename in ['babylon.js', 'babylonjs.loaders.min.js', 'earcut.min.js', 'pep.js']:
    with pkg_resources.resource_stream(
            pkg_resources.Requirement('volmdlr'),
            os.path.join('volmdlr/assets/js/', filename)) as fjs:
        babylon_unpacker_embedded_header += fjs.read().decode('utf-8')

babylon_unpacker_embedded_header += '''
      </script>
</head>
'''


babylon_unpacker_body_template = Template(
'''
<body>
   <canvas id="renderCanvas"></canvas>
   <script type="text/javascript">
      // Get the canvas element from our HTML below
      var canvas = document.querySelector("#renderCanvas");
      // Load the BABYLON 3D engine
      var engine = new BABYLON.Engine(canvas, true);

      var babylon_data = $babylon_data;
      var max_length = babylon_data['max_length'];

      // -------------------------------------------------------------
      // Here begins a function that we will 'call' just after it's built
      var createScene = function () {
        // This creates a basic Babylon Scene object (non-mesh)
        var scene = new BABYLON.Scene(engine);
        scene.useRightHandedSystem = true;
        scene.clearColor = new BABYLON.Color4(.9, .9, .9, .9);
      	var camera = new BABYLON.ArcRotateCamera("ArcRotateCamera",
                                              0, 0, 2*babylon_data['max_length'],
                                              new BABYLON.Vector3(babylon_data['center'][0],
                                                                  babylon_data['center'][1],
                                                                  babylon_data['center'][2]), scene);
      	camera.wheelPrecision=50./babylon_data['max_length']
      	camera.pinchPrecision=50./babylon_data['max_length']
      	camera.panningSensibility=800./babylon_data['max_length'];
      	camera.minZ=0.01*babylon_data['max_length'];
      	camera.attachControl(canvas);
      	camera.inertia = 0;
      	camera.panningInertia = 0;
      	// cam.mode = BABYLON.Camera.ORTHOGRAPHIC_CAMERA;
      	camera.upVector = new BABYLON.Vector3(0, 0, 1);
      	camera.lowerBetaLimit = null;
        camera.upperBetaLimit = null;
        camera.checkCollisions = false;
        camera.lowerRadiusLimit = 0.01*babylon_data['max_length'];
        scene.lastEdgewidthUpdate = Date.now();


        camera.onViewMatrixChangedObservable.add(() => {
            if ((Date.now() - scene.lastEdgewidthUpdate) > 1000){
                scene.lastEdgewidthUpdate = Date.now();
                for (mesh of scene.meshes){
                    var dist = BABYLON.Vector3.Distance(camera.position, mesh.position);
                    mesh.edgesWidth = dist*0.1;
                }
            }
         })

      	var light1 = new BABYLON.HemisphericLight("light1", new BABYLON.Vector3(-1, -1, -1), scene);
      	light1.intensity=0.5;
      	light1.specular = new BABYLON.Color3(0, 0, 0);

      	// var light2 = new BABYLON.SpotLight("Spot0", new BABYLON.Vector3(0, 30, -10), new BABYLON.Vector3(0, -1, 0), 0.8, 2, scene);
      	// light2.diffuse = new BABYLON.Color3(1, 1, 1);
      	// light2.specular = new BABYLON.Color3(1, 1, 1);
        var light2 = new BABYLON.PointLight("light2", new BABYLON.Vector3(0, 0, 0), scene);
        light2.specular = new BABYLON.Color3(0, 0, 0);
        light2.intensity = 0.3;
        light2.parent = camera;

        var light3 = new BABYLON.HemisphericLight("light3", new BABYLON.Vector3(1, 1, 1), scene);
        light3.specular = new BABYLON.Color3(0, 0, 0);
        light3.intensity = 0.50;


        var showAxis = function (size) {
            var makeTextPlane = function (text, color, size) {
                var dynamicTexture = new BABYLON.DynamicTexture("DynamicTexture", 50, scene, true);
                dynamicTexture.hasAlpha = true;
                dynamicTexture.drawText(text, 5, 40, "bold 36px Arial", color, "transparent", true);
                var plane = new BABYLON.Mesh.CreatePlane("TextPlane", size, scene, true);
                plane.material = new BABYLON.StandardMaterial("TextPlaneMaterial", scene);
                plane.material.backFaceCulling = false;
                plane.material.specularColor = new BABYLON.Color3(0, 0, 0);
                plane.material.diffuseTexture = dynamicTexture;
                return plane;
            };
            var axisX = BABYLON.Mesh.CreateLines("axisX", [
                new BABYLON.Vector3.Zero(), new BABYLON.Vector3(size, 0, 0), new BABYLON.Vector3(size * 0.95, 0.05 * size, 0),
                new BABYLON.Vector3(size, 0, 0), new BABYLON.Vector3(size * 0.95, -0.05 * size, 0)
            ], scene);
            axisX.color = new BABYLON.Color3(1, 0, 0);
            var xChar = makeTextPlane("X", "red", size / 10);
            xChar.position = new BABYLON.Vector3(0.9 * size, -0.05 * size, 0);
            var axisY = BABYLON.Mesh.CreateLines("axisY", [
                new BABYLON.Vector3.Zero(), new BABYLON.Vector3(0, size, 0), new BABYLON.Vector3(-0.05 * size, size * 0.95, 0),
                new BABYLON.Vector3(0, size, 0), new BABYLON.Vector3(0.05 * size, size * 0.95, 0)
            ], scene);
            axisY.color = new BABYLON.Color3(0, 1, 0);
            var yChar = makeTextPlane("Y", "green", size / 10);
            yChar.position = new BABYLON.Vector3(0, 0.9 * size, -0.05 * size);
            var axisZ = BABYLON.Mesh.CreateLines("axisZ", [
                new BABYLON.Vector3.Zero(), new BABYLON.Vector3(0, 0, size), new BABYLON.Vector3(0, -0.05 * size, size * 0.95),
                new BABYLON.Vector3(0, 0, size), new BABYLON.Vector3(0, 0.05 * size, size * 0.95)
            ], scene);
            axisZ.color = new BABYLON.Color3(0, 0, 1);
            var zChar = makeTextPlane("Z", "blue", size / 10);
            zChar.position = new BABYLON.Vector3(0, 0.05 * size, 0.9 * size);
        };

        showAxis(1);


        var meshes = []
        for (let mesh_data of babylon_data['meshes']){
          var mesh = new BABYLON.Mesh(mesh_data['name'], scene);
          meshes.push(mesh);

          var normals = [];
          var vertexData = new BABYLON.VertexData();
          BABYLON.VertexData.ComputeNormals(mesh_data['positions'], mesh_data['indices'], normals);

          vertexData.positions = mesh_data['positions'];
          vertexData.indices = mesh_data['indices'];
          vertexData.normals = normals;
          vertexData.applyToMesh(mesh);
          mesh.enableEdgesRendering(0.9);
          mesh.edgesWidth = max_length*0.1;
          mesh.edgesColor = new BABYLON.Color4(0, 0, 0, 0.6);
          var mat = new BABYLON.StandardMaterial("material", scene);
          // mat.diffuseColor = BABYLON.Color3.Green();
          // mat.specularColor = new BABYLON.Color3(0.5, 0.6, 0.87);
          // mat.emissiveColor = new BABYLON.Color3(1, 1, 1);
          // mat.ambientColor = new BABYLON.Color3(0.23, 0.98, 0.53);
          mat.backFaceCulling = false;
          mesh.material = mat;
          mat.diffuseColor = new BABYLON.Color3(mesh_data['color'][0],
                                                mesh_data['color'][1],
                                                mesh_data['color'][2]);
          mat.alpha = mesh_data['alpha'];

        }


        if (babylon_data['steps']){


          var n_primitives = babylon_data['meshes'].length;
          var n_steps = babylon_data['steps'].length;

          var advancedTexture = BABYLON.GUI.AdvancedDynamicTexture.CreateFullscreenUI("UI");
          var animation_stopped = false;

          var buttonHeightInPixels = 50;
          var buttonWidthInPixels = 150;

          // Set container and its size and position
          let heightInPixels = buttonHeightInPixels*7;
          let widthInPixels = buttonWidthInPixels;

          let topInPixels = -canvas.height/2 + heightInPixels/2;
          let leftInPixels = -canvas.width/2 + widthInPixels/2;

          var buttonsContainer = new BABYLON.GUI.StackPanel("buttons_panel");
          buttonsContainer.background = "#263238";
          buttonsContainer.color = "white";
          buttonsContainer.height = ""+heightInPixels+"px";
          buttonsContainer.width = ""+widthInPixels+"px";
          buttonsContainer.top = ""+topInPixels+"px";
          buttonsContainer.left = ""+leftInPixels+"px";

          var step_number_label = new BABYLON.GUI.TextBlock();
          step_number_label.text = "Step n°1/"+n_steps;
          step_number_label.width = ""+buttonWidthInPixels+"px";
          step_number_label.height = ""+buttonHeightInPixels+"px";
          buttonsContainer.addControl(step_number_label);


          var start_stop_button = BABYLON.GUI.Button.CreateSimpleButton("animation", "Stop/Resume animation");
          start_stop_button.width = ""+buttonWidthInPixels+"px";
          start_stop_button.height = ""+buttonHeightInPixels+"px";
          start_stop_button.onPointerUpObservable.add(function(){animation_stopped = !animation_stopped});
          buttonsContainer.addControl(start_stop_button);

          var first_step_button = BABYLON.GUI.Button.CreateSimpleButton("animation", "First step");
          first_step_button.width = ""+buttonWidthInPixels+"px";
          first_step_button.height = ""+buttonHeightInPixels+"px";
          first_step_button.onPointerUpObservable.add(function(){animation_stopped=true; iframe=0; showStep(Math.floor(0))});
          buttonsContainer.addControl(first_step_button);

          var previous_step_button = BABYLON.GUI.Button.CreateSimpleButton("animation", "Previous step");
          previous_step_button.width = ""+buttonWidthInPixels+"px";
          previous_step_button.height = ""+buttonHeightInPixels+"px";
          previous_step_button.onPointerUpObservable.add(function(){animation_stopped=true; iframe-=10; showStep(Math.floor(iframe/10))});
          buttonsContainer.addControl(previous_step_button);

          var next_step_button = BABYLON.GUI.Button.CreateSimpleButton("animation", "Next step");
          next_step_button.width = ""+buttonWidthInPixels+"px";
          next_step_button.height = ""+buttonHeightInPixels+"px";
          next_step_button.onPointerUpObservable.add(function(){animation_stopped=true; iframe+=10; showStep(Math.floor(iframe/10))});
          buttonsContainer.addControl(next_step_button);

          var last_step_button = BABYLON.GUI.Button.CreateSimpleButton("animation", "Last step");
          last_step_button.width = ""+buttonWidthInPixels+"px";
          last_step_button.height = ""+buttonHeightInPixels+"px";
          last_step_button.onPointerUpObservable.add(function(){animation_stopped=true; iframe=10*(n_steps-1); showStep(Math.floor(iframe/10))});
          buttonsContainer.addControl(last_step_button);

          var step_label = new BABYLON.GUI.TextBlock();
          step_label.text = "Label";
          step_label.width = ""+buttonWidthInPixels+"px";
          step_label.height = ""+buttonHeightInPixels+"px";
          buttonsContainer.addControl(step_label);

          advancedTexture.addControl(buttonsContainer);


          var showStep = function (istep){
            // Parts position update
            step_number_label.text = "Step n°"+(istep+1)+"/"+n_steps;
            for(let iprimitive=0; iprimitive<n_primitives; iprimitive++){
                mesh = meshes[iprimitive];
                mesh.position = new BABYLON.Vector3(babylon_data['steps'][istep][iprimitive]['position'][0],
                                                    babylon_data['steps'][istep][iprimitive]['position'][1],
                                                    babylon_data['steps'][istep][iprimitive]['position'][2]);


                mesh.rotation = BABYLON.Vector3.RotationFromAxis(
                  new BABYLON.Vector3(babylon_data['steps'][istep][iprimitive]['orientations'][0][0],
                                      babylon_data['steps'][istep][iprimitive]['orientations'][0][1],
                                      babylon_data['steps'][istep][iprimitive]['orientations'][0][2]),
                  new BABYLON.Vector3(babylon_data['steps'][istep][iprimitive]['orientations'][1][0],
                                      babylon_data['steps'][istep][iprimitive]['orientations'][1][1],
                                      babylon_data['steps'][istep][iprimitive]['orientations'][1][2]),
                  new BABYLON.Vector3(babylon_data['steps'][istep][iprimitive]['orientations'][2][0],
                                      babylon_data['steps'][istep][iprimitive]['orientations'][2][1],
                                      babylon_data['steps'][istep][iprimitive]['orientations'][2][2]));
            }

            if ('label' in babylon_data['steps'][istep]){
              step_label.text = babylon_data['steps'][istep]['label'];
            }
            else{
              step_label.text = "";
            }

          }

        var iframe = 0;

        scene.registerBeforeRender(function () {
                  if (!animation_stopped){
                    if (iframe % 10 == 0){
                      var istep = iframe / 10;

                      showStep(istep);


                    }
                    iframe++;
                    iframe = iframe % (n_steps*10);
                  }

          });
          }


    	return scene;	  };

      var scene = createScene();

      // Register a render loop to repeatedly render the scene
      engine.runRenderLoop(function () {
         scene.render();
      });
      // Watch for browser/canvas resize events
      window.addEventListener("resize", function () {
         engine.resize();
      });

        //scene.debugLayer.show();

   </script>
</body>

</html>

'''
)