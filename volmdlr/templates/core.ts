export class PlotData {
  plot_datas:any;
  context_show:any;
  context_hidden:any;
  minX:number=0;
  maxX:number=0;
  minY:number=0;
  maxY:number=0;
  init_scale:number;
  scale:number;
  last_mouse1X:number;
  last_mouse1Y:number;
  colour_to_plot_data:any;
  select_on_mouse:any;
  select_on_click:any[]=[];
  color_surface_on_mouse:string='lightskyblue';
  color_surface_on_click:string='blue'

  constructor(public data: any,
              public width: number,
              public height: number,
              public coeff_pixel: number) {
    this.width = width;
    this.height = height;
    this.plot_datas = []
    this.select_on_click = []
    this.colour_to_plot_data = {}
    for (var i = 0; i < data.length; i++) {
      var d = data[i]
      if (d['type'] == 'contour'){
        var a = PlotDataContour2D.deserialize(d)
        this.plot_datas.push(a)
        this.minX = Math.min(this.minX, a.minX)
        this.maxX = Math.max(this.maxX, a.maxX)
        this.minY = Math.min(this.minY, a.minY)
        this.maxY = Math.max(this.maxY, a.maxY)
        this.colour_to_plot_data[a.mouse_selection_color] = a
      } else if (d['type'] == 'point') {
        var b = PlotDataPoint2D.deserialize(d);
        this.plot_datas.push(b);
        this.minX = Math.min(this.minX, b.minX)
        this.maxX = Math.max(this.maxX, b.maxX)
        this.minY = Math.min(this.minY, b.minY)
        this.maxY = Math.max(this.maxY, b.maxY)
        this.colour_to_plot_data[b.mouse_selection_color] = b
      } else if (d['type'] == 'plot') {
        var c = PlotDataScatterPlot.deserialize(d)
        this.plot_datas.push(c)
      }
    }
    this.define_canvas()
    this.mouse_interaction()
  }

  define_canvas() {
    var canvas = document.getElementById('canvas');
    canvas.width = this.width;
		canvas.height = this.height;
    this.context_show = canvas.getContext("2d");

    var hiddenCanvas = document.createElement("canvas");
		hiddenCanvas.width = this.width;
		hiddenCanvas.height = this.height;
    this.context_hidden = hiddenCanvas.getContext("2d");
  }

  draw_initial() {
    this.init_scale = Math.min(this.width/(this.coeff_pixel*this.maxX - this.coeff_pixel*this.minX), this.height/(this.coeff_pixel*this.maxY - this.coeff_pixel*this.minY))
    this.scale = this.init_scale
		this.last_mouse1X = (this.width/2 - (this.coeff_pixel*this.maxX - this.coeff_pixel*this.minX)*this.scale/2)/this.scale - this.coeff_pixel*this.minX
		this.last_mouse1Y = (this.height/2 - (this.coeff_pixel*this.maxY - this.coeff_pixel*this.minY)*this.scale/2)/this.scale - this.coeff_pixel*this.minY
		this.draw(false, 0, this.last_mouse1X, this.last_mouse1Y, this.scale)
    this.draw(true, 0, this.last_mouse1X, this.last_mouse1Y, this.scale)
  }

  draw(hidden, show_state, mvx, mvy, scale) {
    if (hidden) {
      var context = this.context_hidden
    } else {
      var context = this.context_show
    }

    context.clearRect(0, 0, this.width, this.height);
    

    for (var i = 0; i < this.plot_datas.length; i++) {
      var d = this.plot_datas[i]
      if (d['type'] == 'contour') {
        context.beginPath();
        if (hidden) {
          context.fillStyle = d.mouse_selection_color;
        } else {
          context.strokeStyle = d.plot_data_states[show_state].color_line
          context.lineWidth = d.plot_data_states[show_state].stroke_width;
          context.fillStyle = 'white'
          if (d.plot_data_states[show_state].hatching != null) {
            context.fillStyle = context.createPattern(d.plot_data_states[show_state].hatching.canvas_hatching,'repeat');
          }
          if (d.plot_data_states[show_state].color_surface != null) {
            context.fillStyle = d.plot_data_states[show_state].color_surface.color;
          }
          if (this.select_on_mouse == d) {
            context.fillStyle = this.color_surface_on_mouse
          }
          for (var j = 0; j < this.select_on_click.length; j++) {
            var z = this.select_on_click[j]
            if (z == d) {
              context.fillStyle = this.color_surface_on_click
            }
          }
        }
        for (var j = 0; j < d.plot_data_primitives.length; j++) {
          var elem = d.plot_data_primitives[j]
          if (j == 0) {var first_elem = true} else {var first_elem = false}
          elem.draw(context, first_elem,  mvx, mvy, scale)
        }
        context.closePath();
        context.fill();
        
      } else if (d['type'] == 'point'){
        if (hidden) {
          context.fillStyle = d.mouse_selection_color;
        } else {
          context.fillStyle = d.plot_data_states[show_state].point_color.color_fill;
          context.strokeStyle = d.plot_data_states[show_state].point_color.color_stroke;

          if (this.select_on_mouse == d) {
            context.fillStyle = this.color_surface_on_mouse
          }

          for (var j = 0; j < this.select_on_click.length; j++) {
            var z = this.select_on_click[j]
            if (z == d) {
              context.fillStyle = this.color_surface_on_click
            }
          }
        }
        context.beginPath()
        d.draw(context, first_elem,  mvx, mvy, scale, this.init_scale)
        context.closePath();
        context.fill();

      } else if (d['type'] == 'plot'){
        context.beginPath()
        d.draw(context, first_elem, mvx, mvy, scale, this.width, this.height)
        context.closePath()
        context.fill()

      } else {
        throw new Error("Invalid type for plotting. For now, only contours and points can be plotted")
      }
      context.stroke();
    }
  }

  mouse_interaction() {
    var isDrawing = false
		var mouse_mouving = false
		var mouse1X = 0
		var mouse1Y = 0
		var mouse2X = 0
		var mouse2Y = 0
		var mouse3X = 0
		var mouse3Y = 0

    var canvas = document.getElementById('canvas')

    canvas.addEventListener('mousedown', e => {
      mouse1X = e.offsetX;
      mouse1Y = e.offsetY;
      isDrawing = true
    })

    canvas.addEventListener('mousemove', e => {
      if (isDrawing === true) {
				mouse_mouving = true
				mouse2X = e.offsetX;
				mouse2Y = e.offsetY;
        this.draw(false, 0, this.last_mouse1X + mouse2X/this.scale - mouse1X/this.scale, this.last_mouse1Y + mouse2Y/this.scale - mouse1Y/this.scale, this.scale)
        this.draw(true, 0, this.last_mouse1X + mouse2X/this.scale - mouse1X/this.scale, this.last_mouse1Y + mouse2Y/this.scale - mouse1Y/this.scale, this.scale)
			}
      else {
				var mouseX = e.offsetX;
				var mouseY = e.offsetY;
				var col = this.context_hidden.getImageData(mouseX, mouseY, 1, 1).data;
				var colKey = 'rgb(' + col[0] + ',' + col[1] + ',' + col[2] + ')';
        this.select_on_mouse = this.colour_to_plot_data[colKey]
        this.draw(false, 0, this.last_mouse1X, this.last_mouse1Y, this.scale)
			}
    })

    canvas.addEventListener('mouseup', e => {
      if (mouse_mouving) {
				this.last_mouse1X = this.last_mouse1X + mouse2X/this.scale - mouse1X/this.scale
				this.last_mouse1Y = this.last_mouse1Y + mouse2Y/this.scale - mouse1Y/this.scale
			}
      else {
  				var col = this.context_hidden.getImageData(mouse1X, mouse1Y, 1, 1).data;
  				var colKey = 'rgb(' + col[0] + ',' + col[1] + ',' + col[2] + ')';
          var click_plot_data = this.colour_to_plot_data[colKey]
          if (this.is_include(click_plot_data, this.select_on_click)) {
            this.select_on_click = this.remove_selection(click_plot_data, this.select_on_click)
          } else {
            this.select_on_click.push(click_plot_data)
          }
          this.draw(false, 0, this.last_mouse1X, this.last_mouse1Y, this.scale)
  			}
      isDrawing = false
			mouse_mouving = false
    })

    canvas.addEventListener('wheel', e => {
      var event = -e.deltaY/100
      this.scale = this.scale + event
      mouse3X = e.offsetX;
      mouse3Y = e.offsetY;
      this.last_mouse1X = this.last_mouse1X - (mouse3X/(this.scale - event) - mouse3X/this.scale)
      this.last_mouse1Y = this.last_mouse1Y - (mouse3Y/(this.scale - event) - mouse3Y/this.scale)
      this.draw(false, 0, this.last_mouse1X, this.last_mouse1Y, this.scale)
      this.draw(true, 0, this.last_mouse1X, this.last_mouse1Y, this.scale)
    })
  }

  remove_selection(val, list){
    var temp = []
    for (var i = 0; i < list.length; i++) {
      var d = list[i]
      if (val != d) {
        temp.push(d)
      }
    }
    return temp
  }

  is_include(val, list){
    var check = false
    for (var i = 0; i < list.length; i++) {
      var d = list[i]
      if (val == d) {
        check = true
      }
    }
    return check
  }

}

export class PlotDataContour2D {
  minX:number=0;
  maxX:number=0;
  minY:number=0;
  maxY:number=0;
  mouse_selection_color:any;

  constructor(public plot_data_primitives:any,
              public plot_data_states:any,
              public type:string,
              public name:string) {
      for (var i = 0; i < this.plot_data_primitives.length; i++) {
        var d = plot_data_primitives[i]
        this.minX = Math.min(this.minX, d.minX)
        this.maxX = Math.max(this.maxX, d.maxX)
        this.minY = Math.min(this.minY, d.minY)
        this.maxY = Math.max(this.maxY, d.maxY)
      }
      this.mouse_selection_color = genColor()

  }

  public static deserialize(serialized) {
      var temp = serialized['plot_data_states']
      var plot_data_states = []
      for (var i = 0; i < temp.length; i++) {
        var d = temp[i]
        plot_data_states.push(PlotDataState.deserialize(d))
      }
      var temp = serialized['plot_data_primitives']
      var plot_data_primitives = []
      for (var i = 0; i < temp.length; i++) {
        var d = temp[i]
        if (d['type'] == 'line') {
          plot_data_primitives.push(PlotDataLine2D.deserialize(d))
        }
        if (d['type'] == 'circle') {
          plot_data_primitives.push(PlotDataCircle2D.deserialize(d))
        }
        if (d['type'] == 'arc') {
          plot_data_primitives.push(PlotDataArc2D.deserialize(d))
        }

      }
      return new PlotDataContour2D(plot_data_primitives,
                                   plot_data_states,
                                   serialized['type'],
                                   serialized['name']);
  }
}

export class PlotDataLine2D {
  minX:number=0;
  maxX:number=0;
  minY:number=0;
  maxY:number=0;

  constructor(public data:any,
              public plot_data_states:any,
              public type:string,
              public name:string) {
      this.minX = Math.min(this.data[0], this.data[2]);
      this.maxX = Math.max(this.data[0], this.data[2]);
      this.minY = Math.min(this.data[1], this.data[3]);
      this.maxY = Math.max(this.data[1], this.data[3]);
  }

  public static deserialize(serialized) {
      var temp = serialized['plot_data_states']
      var plot_data_states = []
      for (var i = 0; i < temp.length; i++) {
        var d = temp[i]
        plot_data_states.push(PlotDataState.deserialize(d))
      }
      return new PlotDataLine2D(serialized['data'],
                               plot_data_states,
                               serialized['type'],
                               serialized['name']);
  }

  draw(context, first_elem, mvx, mvy, scale) {
    if (first_elem) {
      context.moveTo(scale*(1000*this.data[0]+ mvx), scale*(1000*this.data[1]+ mvy));
    }
    context.lineTo(scale*(1000*this.data[2]+ mvx), scale*(1000*this.data[3]+ mvy));
  }
}

export class PlotDataCircle2D {
  minX:number=0;
  maxX:number=0;
  minY:number=0;
  maxY:number=0;

  constructor(public data:any,
              public cx:number,
              public cy:number,
              public r:number,
              public plot_data_states:PlotDataState[],
              public type:string,
              public name:string) {
      this.minX = this.cx - this.r;
      this.maxX = this.cx + this.r;
      this.minY = this.cy - this.r;
      this.maxY = this.cy + this.r;
              }

  public static deserialize(serialized) {
      var temp = serialized['plot_data_states']
      var plot_data_states = []
      for (var i = 0; i < temp.length; i++) {
        var d = temp[i]
        plot_data_states.push(PlotDataState.deserialize(d))
      }
      return new PlotDataCircle2D(serialized['data'],
                                  serialized['cx'],
                                  serialized['cy'],
                                  serialized['r'],
                                  plot_data_states,
                                  serialized['type'],
                                  serialized['name']);
  }

  draw(context, first_elem, mvx, mvy, scale) {
    context.arc(scale*(1000*this.cx+ mvx), scale*(1000*this.cy+ mvy), scale*1000*this.r, 0, 2*Math.PI);
  }

}

export class PlotDataPoint2D {
  minX:number=0;
  maxX:number=0;
  minY:number=0;
  maxY:number=0;
  mouse_selection_color:any;

  constructor(public data:any,
              public cx:number,
              public cy:number,
              public size:number,
              public plot_data_states:PlotDataState[],
              public type:string,
              public name:string) {
    
    for (var i=0; i<this.plot_data_states.length; i++) {
      var plot = this.plot_data_states[i]
      var point_size = plot.point_size.size
      if (point_size==1||point_size==2||point_size==3||point_size==4) {
        var height = plot.window_size.height
      var width = plot.window_size.width
      } else {
        throw new Error('Invalid point_size')
      }
    }
    this.size = point_size * Math.min(height,width)/200
    this.minX = this.cx - this.size;
    this.maxX = this.cx + this.size;
    this.minY = this.cy - this.size;
    this.maxY = this.cy + this.size;
    this.mouse_selection_color = genColor()
    }

    public static deserialize(serialized) {
      var temp = serialized['plot_data_states']
      var plot_data_states = []
      for (var i = 0; i < temp.length; i++) {
        var d = temp[i]
        plot_data_states.push(PlotDataState.deserialize(d))
      }
      return new PlotDataPoint2D(serialized['data'],
                                  serialized['cx'],
                                  serialized['cy'],
                                  serialized['r'],
                                  plot_data_states,
                                  serialized['type'],
                                  serialized['name']);
    }

    draw(context, first_elem, mvx, mvy, scale, init_scale) {

        for (var i=0; i<this.plot_data_states.length; i++) {
          var shape = this.plot_data_states[i].shape_set.shape
          if (shape == 'circle') {
            context.arc(scale*(1000*this.cx+ mvx), scale*(1000*this.cy+ mvy), init_scale*1000*this.size, 0, 2*Math.PI);
          } else if (shape == 'square') {
            context.rect(scale*(1000*(this.cx - this.size) + mvx),scale*(1000*(this.cy - this.size) + mvy),init_scale*1000*this.size*2, init_scale*1000*this.size*2)
            context.stroke()
          } else if (shape == 'crux') {
            context.rect(scale*(1000*this.cx + mvx), scale*(1000*this.cy + mvy),init_scale*1000*this.size, init_scale*10*this.size)
            context.rect(scale*(1000*this.cx + mvx), scale*(1000*this.cy + mvy),-init_scale*1000*this.size, init_scale*10*this.size)
            context.rect(scale*(1000*this.cx + mvx), scale*(1000*this.cy + mvy),init_scale*10*this.size, init_scale*1000*this.size)
            context.rect(scale*(1000*this.cx + mvx), scale*(1000*this.cy + mvy),init_scale*10*this.size, -init_scale*1000*this.size)
            context.stroke()
          } else {
            throw new Error('Invalid shape for point')
          }
        }
    }
}

export class PlotDataScatterPlot {
  colorStroke:any;
  nb_points_x:number;
  nb_points_y:number;
  constructor(public x_start:number,
                     public x_end:number,
                     public x_step:number,
                     public y_start:number,
                     public y_end:number,
                     public y_step:number, 
                     public name:string, 
                     public type:string, 
                     public plot_data_states:PlotDataState[]) {
    this.nb_points_x = (this.x_end - this.x_start)/this.x_step
    if (this.nb_points_x != Math.abs(Math.floor(this.nb_points_x))) {
      throw new Error('ScatterPlot x_step is not valid')
    }
    this.nb_points_y = (this.y_end - this.y_start)/this.y_step
    if (this.nb_points_y != Math.abs(Math.floor(this.nb_points_y))) {
      throw new Error('ScatterPlot y_step is not valid')
    }

    for (var i=0; i<this.plot_data_states.length; i++) {
      var plot = this.plot_data_states[i]
      this.colorStroke = plot.color_line
    }
  }

  public static deserialize(serialized) {
    var temp = serialized['plot_data_states']
    var plot_data_states = []
    for (var i = 0; i < temp.length; i++) {
      var d = temp[i]
      plot_data_states.push(PlotDataState.deserialize(d))
    }
    return new PlotDataScatterPlot(serialized['x_start'],
                           serialized['x_end'],
                           serialized['x_step'],
                           serialized['y_start'],
                           serialized['y_end'],
                           serialized['y_step'],
                           serialized['name'],
                           serialized['type'],
                           serialized['plot_data_states'])
  }

  draw(context, first_elem, mvx, mvy, scale, width, height) {
    // Dessin du repère
    context.strokeStyle = this.colorStroke
    //Flèches

    context.moveTo(0, 20)
    context.lineTo(10, 0)
    context.moveTo(10, 0)
    context.lineTo(20, 20)
    
    context.moveTo(width - 20, height - 20)
    context.lineTo(width, height - 10)
    context.moveTo(width, height - 10)
    context.lineTo(width - 20, height)

    //Axes
    context.moveTo(10, 0)
    context.lineTo(10, height)

    context.moveTo(0, height - 10)
    context.lineTo(width, height - 10)

    //Graduations
      //pour l'axe des x
    var display_nb_x = 0
    if (this.x_step.toString().split(".").length == 2) {
      display_nb_x = this.x_step.toString().split(".")[1].length;
    }
    for (var i=0; i<this.nb_points_x; i++) {
      context.moveTo(scale*(1000*(this.x_start + i*this.x_step) + mvx), height - 13)
      context.lineTo(scale*(1000*(this.x_start + i*this.x_step) + mvx), height - 7)

      context.font = '12px Arial';
      context.textAlign = 'start';
      context.fillText(Math.round((this.x_start + i*this.x_step)*Math.pow(10,display_nb_x))/Math.pow(10,display_nb_x), scale*(1000*(this.x_start + i*this.x_step) - 20 + mvx), height - 15 )
      
    }
      //pour l'axe des y
    var display_nb_y = 0
    if (this.y_step.toString().split(".").length == 2) {
      display_nb_y = this.x_step.toString().split(".")[1].length;
    }
    for (var i=0; i<this.nb_points_y; i++) {
      context.moveTo(7, scale*(-1000*(this.y_start + i*this.y_step) + mvy))
      context.lineTo(13, scale*(-1000*(this.y_start + i*this.y_step) + mvy))
      context.fillText(Math.round((this.y_start + i*this.y_step)*Math.pow(10,display_nb_y))/Math.pow(10,display_nb_y), 15, scale*(-1000*(this.y_start + i*this.y_step) + mvy))
    }

    context.stroke()
  }
}

export class PlotDataArc2D {
  minX:number=0;
  maxX:number=0;
  minY:number=0;
  maxY:number=0;

  constructor(public cx:number,
              public cy:number,
              public r:number,
              public data:any,
              public angle1:number,
              public angle2:number,
              public plot_data_states:PlotDataState[],
              public type:string,
              public name:string) {
      if((this.cx - this.r) < this.minX){
        this.minX = this.cx - this.r;
      }
      if((this.cx - this.r) > this.maxX){
        this.maxX = this.cx + this.r;
      }
      if((this.cy - this.r) < this.minY){
        this.minY = this.cy - this.r;
      }
      if((this.cy + this.r) > this.maxY){
        this.maxY = this.cy + this.r;
      }
  }

  public static deserialize(serialized) {
      var temp = serialized['plot_data_states']
      var plot_data_states = []
      for (var i = 0; i < temp.length; i++) {
        var d = temp[i]
        plot_data_states.push(PlotDataState.deserialize(d))
      }
      return new PlotDataArc2D(serialized['cx'],
                                  serialized['cy'],
                                  serialized['r'],
                                  serialized['data'],
                                  serialized['angle1'],
                                  serialized['angle2'],
                                  plot_data_states,
                                  serialized['type'],
                                  serialized['name']);
  }

  draw(context, first_elem, mvx, mvy, scale) {
    var ptsa = []
    for (var l = 0; l < this.data.length; l++) {
      ptsa.push(scale*(1000*this.data[l]['x']+ mvx))
      ptsa.push(scale*(1000*this.data[l]['y']+ mvy))
    }
    var tension = 0.4
    var isClosed = false
    var numOfSegments = 16
    drawLines(context, getCurvePoints(ptsa, tension, isClosed, numOfSegments));
  }
}

export class PlotDataState {

  constructor(public color_surface:ColorSurfaceSet,
              public color_map:any,
              public hatching:HatchingSet,
              public opacity:number,
              public dash:any,
              public marker:any,
              public color_line:any,
              public shape_set:PointShapeSet,
              public point_size:PointSizeSet,
              public point_color:PointColorSet,
              public window_size:WindowSizeSet,
              public stroke_width:any,
              public name:any,) {}

  public static deserialize(serialized) {
      var color_surface = null
      if (serialized['color_surface'] != null) {
        color_surface = ColorSurfaceSet.deserialize(serialized['color_surface'])
      }
      var hatching = null
      if (serialized['hatching'] != null) {
        hatching = HatchingSet.deserialize(serialized['hatching'])
      }
      var shape_set = null
      if (serialized['shape_set'] != null) {
        shape_set = PointShapeSet.deserialize(serialized['shape_set'])
      }
      var window_size = null
      if(serialized['window_size'] != null) {
        window_size = WindowSizeSet.deserialize(serialized['window_size'])
      }
      var point_size = null
      if (serialized['point_size'] != null) {
        point_size = PointSizeSet.deserialize(serialized['point_size'])
      }
      var point_color = null
      if (serialized['point_color'] != null) {
        point_color = PointColorSet.deserialize(serialized['point_color'])
      }
      return new PlotDataState(color_surface,
                               serialized['color_map'],
                               hatching,
                               serialized['opacity'],
                               serialized['dash'],
                               serialized['marker'],
                               serialized['color_line'],
                               shape_set,
                               point_size,
                               point_color,
                               window_size,
                               serialized['stroke_width'],
                               serialized['name']);
  }
}

export class ColorSurfaceSet {

  constructor(public name:string,
              public color:any) {}

  public static deserialize(serialized) {
      return new ColorSurfaceSet(serialized['name'],
                               serialized['color']);
  }
}

export class PointShapeSet {
  constructor(public name:string, public shape:any){}

  public static deserialize(serialized) {
    return new PointShapeSet(serialized['name'],
                             serialized['shape']);
  }
}

export class PointSizeSet {
  constructor(public name:string, public size:number) {}

  public static deserialize(serialized) {
    return new PointSizeSet(serialized['name'],
                            serialized['size'])
  }
}

export class PointColorSet {
  constructor(public name:string, public color_fill:string, public color_stroke:string) {}

  public static deserialize(serialized) {
    return new PointColorSet(serialized['name'],
                             serialized['color_fill'],
                             serialized['color_stroke'])
  }
}

export class WindowSizeSet {
  constructor(public name:string, public height:number,public width:number){
  }

  public static deserialize(serialized) {
    return new WindowSizeSet(serialized['name'],
                             serialized['height'],
                             serialized['width']);
  }
}

export class HatchingSet {
  canvas_hatching:any;

  constructor(public name:string,
              public stroke_width:number,
              public hatch_spacing:number) {
      this.canvas_hatching = this.generate_canvas()
  }

  public static deserialize(serialized) {
      return new HatchingSet(serialized['name'],
                             serialized['stroke_width'],
                             serialized['hatch_spacing']);
  }

  generate_canvas() {
    var nb_hatch = 20
    var max_size = nb_hatch*this.hatch_spacing

    var p_hatch = document.createElement("canvas");
    p_hatch.width = max_size;
    p_hatch.height = max_size;
    var pctx = p_hatch.getContext("2d");
    pctx.lineCap = 'square';
    pctx.strokeStyle = 'black'
    pctx.lineWidth = this.stroke_width;
    pctx.beginPath();
    var pos_x = - Math.pow(Math.pow(max_size,2)/2, 0.5)
    var pos_y = Math.pow(Math.pow(max_size,2)/2, 0.5)
    for (var i = 0; i <= 2*nb_hatch; i++) {
      pos_x = pos_x + this.hatch_spacing
      pos_y = pos_y - this.hatch_spacing
      pctx.moveTo(pos_x, pos_y);
      pctx.lineTo(pos_x + max_size, pos_y + max_size);
    }
    pctx.stroke();
    return p_hatch
  }
}


function drawLines(ctx, pts) {
    // ctx.moveTo(pts[0], pts[1]);
    for(var i=2; i<pts.length-1; i+=2) ctx.lineTo(pts[i], pts[i+1]);
}

function getCurvePoints(pts, tension, isClosed, numOfSegments) {

    // use input value if provided, or use a default value
    tension = (typeof tension != 'undefined') ? tension : 0.5;
    isClosed = isClosed ? isClosed : false;
    numOfSegments = numOfSegments ? numOfSegments : 16;

    var _pts = [], res = [],    // clone array
        x, y,           // our x,y coords
        t1x, t2x, t1y, t2y, // tension vectors
        c1, c2, c3, c4,     // cardinal points
        st, t, i;       // steps based on num. of segments

    // clone array so we don't change the original
    //
    _pts = pts.slice(0);

    // The algorithm require a previous and next point to the actual point array.
    // Check if we will draw closed or open curve.
    // If closed, copy end points to beginning and first points to end
    // If open, duplicate first points to befinning, end points to end
    if (isClosed) {
        _pts.unshift(pts[pts.length - 1]);
        _pts.unshift(pts[pts.length - 2]);
        _pts.unshift(pts[pts.length - 1]);
        _pts.unshift(pts[pts.length - 2]);
        _pts.push(pts[0]);
        _pts.push(pts[1]);
    }
    else {
        _pts.unshift(pts[1]);   //copy 1. point and insert at beginning
        _pts.unshift(pts[0]);
        _pts.push(pts[pts.length - 2]); //copy last point and append
        _pts.push(pts[pts.length - 1]);
    }

    // ok, lets start..

    // 1. loop goes through point array
    // 2. loop goes through each segment between the 2 pts + 1e point before and after
    for (i=2; i < (_pts.length - 4); i+=2) {
        for (t=0; t <= numOfSegments; t++) {

            // calc tension vectors
            t1x = (_pts[i+2] - _pts[i-2]) * tension;
            t2x = (_pts[i+4] - _pts[i]) * tension;

            t1y = (_pts[i+3] - _pts[i-1]) * tension;
            t2y = (_pts[i+5] - _pts[i+1]) * tension;

            // calc step
            st = t / numOfSegments;

            // calc cardinals
            c1 =   2 * Math.pow(st, 3)  - 3 * Math.pow(st, 2) + 1;
            c2 = -(2 * Math.pow(st, 3)) + 3 * Math.pow(st, 2);
            c3 =       Math.pow(st, 3)  - 2 * Math.pow(st, 2) + st;
            c4 =       Math.pow(st, 3)  -     Math.pow(st, 2);

            // calc x and y cords with common control vectors
            x = c1 * _pts[i]    + c2 * _pts[i+2] + c3 * t1x + c4 * t2x;
            y = c1 * _pts[i+1]  + c2 * _pts[i+3] + c3 * t1y + c4 * t2y;

            //store points in array
            res.push(x);
            res.push(y);

        }
    }

    return res;
}

var nextCol = 1;
function genColor(){
  var ret = [];
  // via http://stackoverflow.com/a/15804183
  if(nextCol < 16777215){
    ret.push(nextCol & 0xff); // R
    ret.push((nextCol & 0xff00) >> 8); // G
    ret.push((nextCol & 0xff0000) >> 16); // B

    nextCol += 1;
  }
  var col = "rgb(" + ret.join(',') + ")";
  return col;
}