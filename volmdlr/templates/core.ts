export class PlotData {
  context_show:any;
  context_hidden:any;
  minX:number=0;
  maxX:number=0;
  minY:number=0;
  maxY:number=0;
  init_scale:number;
  scale:number;
  scaleX:number;
  scaleY:number;
  scroll_x:number=0;
  scroll_y:number=0;
  last_mouse1X:number;
  last_mouse1Y:number;
  colour_to_plot_data:any={};
  select_on_mouse:any;
  select_on_click:any[]=[];
  tooltip_list:any[]=[];
  color_surface_on_mouse:string='lightskyblue';
  color_surface_on_click:string='blue';
  zoom_rect_ON:boolean=true;
  zoom_rect_x:number;
  zoom_rect_y:number;
  zoom_rect_w:number;
  zoom_rect_h:number;
  zw_ON:boolean=true;
  zw_bool:boolean;
  zw_x:number;
  zw_y:number;
  zw_w:number;
  zw_h:number;
  reset_ON:boolean=true;
  reset_rect_x:number;
  reset_rect_y:number;
  reset_rect_w:number;
  reset_rect_h:number;
  select_bool:boolean;
  select_x:number;
  select_y:number;
  select_w:number;
  select_h:number;

  define_canvas(width, height) {
    var canvas = document.getElementById('canvas');
    canvas.width = width;
		canvas.height = height;
    this.context_show = canvas.getContext("2d");

    var hiddenCanvas = document.createElement("canvas");
		hiddenCanvas.width = width;
		hiddenCanvas.height = height;
    this.context_hidden = hiddenCanvas.getContext("2d");
  }

  draw_initial(width, height, coeff_pixel) {
    this.init_scale = Math.min(width/(coeff_pixel*this.maxX - coeff_pixel*this.minX), height/(coeff_pixel*this.maxY - coeff_pixel*this.minY));
    this.scale = this.init_scale;
    this.scaleX = this.init_scale;
    this.scaleY = this.init_scale;
		this.last_mouse1X = (width/2 - (coeff_pixel*this.maxX - coeff_pixel*this.minX)*this.scale/2)/this.scale - coeff_pixel*this.minX;
    this.last_mouse1Y = (height/2 - (coeff_pixel*this.maxY - coeff_pixel*this.minY)*this.scale/2)/this.scale - coeff_pixel*this.minY;
  }

  draw(hidden, show_state, mvx, mvy, scaleX, scaleY, width, height, plot_datas) {
    if (hidden) {
      var context = this.context_hidden;
    } else {
      var context = this.context_show;
    }
    context.clearRect(0, 0, width, height);

    for (var i = 0; i < plot_datas.length; i++) {
      var d = plot_datas[i];

      if (d['type'] == 'contour') {
        context.beginPath();
        if (hidden) {
          context.fillStyle = d.mouse_selection_color;
        } else {
          context.strokeStyle = d.plot_data_states[show_state].color_line;
          context.lineWidth = d.plot_data_states[show_state].stroke_width;
          context.fillStyle = 'white';
          if (d.plot_data_states[show_state].hatching != null) {
            context.fillStyle = context.createPattern(d.plot_data_states[show_state].hatching.canvas_hatching,'repeat');
          }
          if (d.plot_data_states[show_state].color_surface != null) {
            context.fillStyle = d.plot_data_states[show_state].color_surface.color;
          }
          if (this.select_on_mouse == d) {
            context.fillStyle = this.color_surface_on_mouse;
          }
          for (var j = 0; j < this.select_on_click.length; j++) {
            var z = this.select_on_click[j];
            if (z == d) {
              context.fillStyle = this.color_surface_on_click;
            }
          }
        }
        for (var j = 0; j < d.plot_data_primitives.length; j++) {
          var elem = d.plot_data_primitives[j];
          if (j == 0) {var first_elem = true} else {var first_elem = false}
          elem.draw(context, first_elem,  mvx, mvy, scaleX, scaleY);
        }
        context.fill();
        context.stroke();
        context.closePath();
      
      } else if (d['type'] == 'point') {

        if (hidden) {
          context.fillStyle = d.mouse_selection_color;
        } else {
          context.fillStyle = d.plot_data_states[show_state].point_color.color_fill;
          context.lineWidth = d.plot_data_states[show_state].stroke_width;
          context.strokeStyle = d.plot_data_states[show_state].point_color.color_stroke;

          if (this.select_on_mouse == d) {
            context.fillStyle = this.color_surface_on_mouse;
          }
          for (var j = 0; j < this.select_on_click.length; j++) {
            var z = this.select_on_click[j];
            var shape = d.plot_data_states[show_state].shape_set.shape;

            if (z == d) {
              if (shape == 'crux') {
                context.strokeStyle = this.color_surface_on_click;
              } else {
                context.fillStyle = this.color_surface_on_click    ;         
              }
            }
          }
        }
        var x = scaleX*(1000*d.cx+ mvx);
        var y = scaleY*(1000*d.cy + mvy);
        var length = 1000*d.size*this.init_scale;

        var is_inside_canvas = ((x + length>=0) && (x - length <= width) && (y + length >= 0) && (y - length <= height));

        if (is_inside_canvas === true) {
          context.lineWidth = 1
          context.beginPath();
          d.draw(context, this.context_hidden, mvx, mvy, scaleX, scaleY, this.init_scale);
          context.fill();
          context.closePath();
        }

      } else if (d['type'] == 'plot'){
        context.beginPath();
        d.draw(context, mvx, mvy, scaleX, scaleY, width, height, this.init_scale, this.minX, this.maxX, this.minY, this.maxY, this.scroll_x, this.scroll_y);
        context.closePath();
        context.fill();

      }
    }
      
      //Drawing the tooltips
      for (var i=0; i<this.tooltip_list.length; i++) {
        if (!(typeof this.tooltip_list[i] === "undefined")) {
          this.tooltip(context, this.tooltip_list[i], scaleX, scaleY, mvx, mvy,width, height);
        }
      }
  
      //Drawing the zooming button
      if (this.zoom_rect_ON) {
        this.zoom_rect_x = width - 45;
        this.zoom_rect_y = 10;
        this.zoom_rect_w = 35;
        this.zoom_rect_h = 25;
        this.zoom_button(context, this.zoom_rect_x, this.zoom_rect_y, this.zoom_rect_w, this.zoom_rect_h, width, height);
      }
      
      //Drawing the button for zooming window selection
      if (this.zw_ON) {
        this.zw_x = width - 45;
        this.zw_y = 70;
        this.zw_w = 35;
        this.zw_h = 30;
        this.zoom_window_button(context, this.zw_x,this.zw_y,this.zw_w,this.zw_h, width, height);
      }
  
      //Drawing the reset button
      if (this.reset_ON) {
        this.reset_rect_x = width - 45;
        this.reset_rect_y = 110;
        this.reset_rect_w = 35;
        this.reset_rect_h = 30;
        this.reset_button(context, this.reset_rect_x, this.reset_rect_y, this.reset_rect_w, this.reset_rect_h, width, height);
      }
  
      this.select_x = width - 45;
      this.select_y = 150;
      this.select_w = 35;
      this.select_h = 30;
      this.selection_button(context, this.select_x, this.select_y, this.select_w, this.select_h, width, height);
    
  }

  tooltip(context, point, scaleX, scaleY, mvx, mvy, width, height) {
    context.beginPath();
    var cx = point.cx;
    var cy = point.cy;

    var rect_w = this.init_scale*300;
    var rect_h = this.init_scale*200;
    var rect_radius = this.init_scale*40;
    var rect_x = scaleX*(1000*cx + mvx) + this.init_scale*80;
    var rect_y = scaleY*(1000*cy + mvy) - 1/2*rect_h;

    if (rect_x + rect_w  > width) {
      rect_x = scaleX*(1000*cx + mvx) - this.init_scale*80 - rect_w;
    }
    if (rect_y < 0) {
      rect_y = scaleY*(1000*cy + mvy);
    }
    if (rect_y + rect_h > height) {
      rect_y = scaleY*(1000*cy + mvy) - rect_h;
    }

    Shape.roundRect(rect_x, rect_y, rect_w, rect_h, rect_radius, context)
    context.strokeStyle = 'black';
    context.fillStyle = 'lightblue';
    context.stroke();
    context.fill();

    var coordinate_size = this.init_scale*50;
    context.font = coordinate_size.toString() + 'px Arial';
    context.fillStyle = 'black';
    context.textAlign = 'center';
    var int_nb_digit_x = cx.toString().split('.')[0].length;
    var nb_digits_x = Math.max(0, 5 - int_nb_digit_x);
    var int_nb_digit_y = (-cy).toString().split('.')[0].length;
    var nb_digits_y = Math.max(0, 5 - int_nb_digit_y);
    var round_cx = MyMath.round(cx,nb_digits_x);
    var round_cy = MyMath.round(cy,nb_digits_y);

    var x_middle = rect_x + 1/2*rect_w;
    var y_middle = rect_y + 1/2*rect_h + rect_radius;
    context.fillText('x = ' + round_cx.toString(), x_middle, y_middle - this.init_scale*65);
    context.fillText('y = ' + (-round_cy).toString(), x_middle, y_middle + this.init_scale*15);
    context.closePath();
  }

  zoom_button(context, x, y, w, h,width, height) {
    if ((x<0) || (x+h>width) || (y<0) || (y+2*h>height)) {
      throw new Error("Invalid x or y, the zoom button is out of the canvas");
    }
    context.beginPath();
    context.lineWidth = "2";
    context.fillStyle = 'white';
    context.rect(x, y, w, h);
    context.rect(x, y+h, w, h);
    context.moveTo(x, y+h);
    context.lineTo(x+w, y+h);
    Shape.crux(context, x+w/2, y+h/2, h/3);
    context.moveTo(x + w/2 - h/3, y + 3*h/2);
    context.lineTo(x + w/2 + h/3, y + 3*h/2);
    context.fill();
    context.stroke();
    context.closePath();
  }

  zoom_window_button(context, x, y, w, h, width, height) {
    if ((x<0) || (x+h>width) || (y<0) || (y+h>height)) {
      throw new Error("Invalid x or y, the zoom window button is out of the canvas");
    }
    if (this.zw_bool) {
      Shape.createButton(x, y, w, h, context, "Z ON");
    } else {
      Shape.createButton(x, y, w, h, context, "Z OFF");
    }
    
  }

  reset_button(context, x, y, w, h, width, height) {
    if ((x<0) || (x+h>width) || (y<0) || (y+h>height)) {
      throw new Error("Invalid x or y, the reset button is out of the canvas");
    }
    Shape.createButton(x, y, w, h, context, "Reset")
  }

  selection_button(context, x, y, w, h, width, height) {
    if ((x<0) || (x+h>width) || (y<0) || (y+h>height)) {
      throw new Error("Invalid x or y, the selection button is out of the canvas");
    }
    if (this.select_bool) {
      Shape.createButton(x, y, w, h, context, "S ON")
    } else {
      Shape.createButton(x, y, w, h, context, "S OFF")
    }
  }

  mouse_interaction(width, height, plot_datas, coeff_pixel) {
    var isDrawing = false;
    var mouse_mouving = false;
    var mouse1X = 0;
    var mouse1Y = 0;
    var mouse2X = 0;
    var mouse2Y = 0;
    var mouse3X = 0;
    var mouse3Y = 0;

    var canvas = document.getElementById('canvas');

    canvas.addEventListener('mousedown', e => {
      mouse1X = e.offsetX;
      mouse1Y = e.offsetY;
      mouse2X = e.offsetX;
      mouse2Y = e.offsetY;
      isDrawing = true;
    })

    canvas.addEventListener('mousemove', e => {
      if ((isDrawing === true) && !(this.zw_bool||this.select_bool)) {
        mouse_mouving = true;
        mouse2X = e.offsetX;
        mouse2Y = e.offsetY;
        this.draw(false, 0, this.last_mouse1X + mouse2X/this.scaleX - mouse1X/this.scaleX, this.last_mouse1Y + mouse2Y/this.scaleY - mouse1Y/this.scaleY, this.scaleX, this.scaleY, width, height, plot_datas);
        this.draw(true, 0, this.last_mouse1X + mouse2X/this.scaleX - mouse1X/this.scaleX, this.last_mouse1Y + mouse2Y/this.scaleY - mouse1Y/this.scaleY, this.scaleX, this.scaleY, width, height, plot_datas);
      
      } else if ((isDrawing === true) && (this.zw_bool||this.select_bool)) {
        mouse_mouving = true;
        mouse2X = e.offsetX;
        mouse2Y = e.offsetY;
        this.draw(false, 0, this.last_mouse1X, this.last_mouse1Y, this.scaleX, this.scaleY, width, height, plot_datas);
        this.draw(true, 0, this.last_mouse1X, this.last_mouse1Y, this.scaleX, this.scaleY, width, height, plot_datas);
        this.context_show.beginPath();
        this.context_show.rect(mouse1X, mouse1Y, mouse2X - mouse1X, mouse2Y - mouse1Y);
        this.context_show.stroke();
        this.context_show.closePath();
        this.context_hidden.beginPath();
        this.context_hidden.rect(mouse1X, mouse1Y, mouse2X - mouse1X, mouse2Y - mouse1Y);
        this.context_hidden.stroke();
        this.context_hidden.closePath(); 
      } else {
        var mouseX = e.offsetX;
        var mouseY = e.offsetY;
        var col = this.context_hidden.getImageData(mouseX, mouseY, 1, 1).data;
        var colKey = 'rgb(' + col[0] + ',' + col[1] + ',' + col[2] + ')';
        this.select_on_mouse = this.colour_to_plot_data[colKey];
        this.draw(false, 0, this.last_mouse1X, this.last_mouse1Y, this.scaleX, this.scaleY, width, height, plot_datas);
      }
    })

    canvas.addEventListener('mouseup', e => {

      var click_on_plus = Shape.Is_in_rect(mouse1X, mouse1Y, this.zoom_rect_x, this.zoom_rect_y, this.zoom_rect_w, this.zoom_rect_h);
      var click_on_minus = Shape.Is_in_rect(mouse1X, mouse1Y, this.zoom_rect_x, this.zoom_rect_y + this.zoom_rect_h, this.zoom_rect_w, this.zoom_rect_h);
      var click_on_zoom_window = Shape.Is_in_rect(mouse1X, mouse1Y, this.zw_x, this.zw_y, this.zw_w, this.zw_h);
      var click_on_reset = Shape.Is_in_rect(mouse1X, mouse1Y, this.reset_rect_x, this.reset_rect_y, this.reset_rect_w, this.reset_rect_h);
      var is_rect_big_enough = (Math.abs(mouse2X - mouse1X)>40) && (Math.abs(mouse2Y - mouse1Y)>30);
      var click_on_select = Shape.Is_in_rect(mouse1X, mouse1Y, this.select_x, this.select_y, this.select_w, this.select_h);
      var click_on_button = click_on_plus || click_on_minus || click_on_zoom_window || click_on_reset || click_on_select;

      if (mouse_mouving) {
          if (this.zw_bool && is_rect_big_enough) {
            var zoom_coeff_x = width/Math.abs(mouse2X - mouse1X);
            var zoom_coeff_y = height/Math.abs(mouse2Y - mouse1Y);
            this.last_mouse1X = this.last_mouse1X - Math.min(mouse1X, mouse2X)/this.scaleX
            this.last_mouse1Y = this.last_mouse1Y - Math.min(mouse1Y,mouse2Y)/this.scaleY
            this.scaleX = this.scaleX*zoom_coeff_x;
            this.scaleY = this.scaleY*zoom_coeff_y;
            this.draw(false, 0, this.last_mouse1X, this.last_mouse1Y, this.scaleX, this.scaleY, width, height, plot_datas);
            this.draw(true, 0, this.last_mouse1X, this.last_mouse1Y, this.scaleX, this.scaleY, width, height, plot_datas);

          } else if (this.select_bool) {
            for (var i=0; i<plot_datas.length; i++) {
              var d = plot_datas[i];
              var in_rect = Shape.Is_in_rect(this.scaleX*(1000*d.cx + this.last_mouse1X),this.scaleY*(1000*d.cy + this.last_mouse1Y), Math.min(mouse1X, mouse2X), Math.min(mouse1Y, mouse2Y), Math.abs(mouse2X - mouse1X), Math.abs(mouse2Y - mouse1Y));
              if ((d['type']=="point") && (in_rect === true) && !(this.is_include(d, this.select_on_click))) {
                this.select_on_click.push(d);
              }
            }
            this.draw(false, 0, this.last_mouse1X, this.last_mouse1Y, this.scaleX, this.scaleY, width, height, plot_datas);
            this.draw(true, 0, this.last_mouse1X, this.last_mouse1Y, this.scaleX, this.scaleY, width, height, plot_datas);

          } else {
            this.last_mouse1X = this.last_mouse1X + mouse2X/this.scaleX - mouse1X/this.scaleX;
            this.last_mouse1Y = this.last_mouse1Y + mouse2Y/this.scaleY - mouse1Y/this.scaleY;
          }

      } else {
          var col = this.context_hidden.getImageData(mouse1X, mouse1Y, 1, 1).data;
          var colKey = 'rgb(' + col[0] + ',' + col[1] + ',' + col[2] + ')';
          var click_plot_data = this.colour_to_plot_data[colKey];
          if (this.is_include(click_plot_data, this.select_on_click)) {
            this.select_on_click = this.remove_selection(click_plot_data, this.select_on_click);
          } else {
            this.select_on_click.push(click_plot_data);
          }
          if (this.is_include(click_plot_data, this.tooltip_list) && (!this.is_include(click_plot_data, this.select_on_click))) {
            this.tooltip_list = this.remove_selection(click_plot_data, this.tooltip_list);
          } else if (!this.is_include(click_plot_data, this.tooltip_list) && this.is_include(click_plot_data, this.select_on_click)){
            this.tooltip_list.push(click_plot_data);
          }

          

          if (this.contains_undefined(this.select_on_click) && !click_on_button) {
            this.select_on_click = [];
            this.tooltip_list = [];
          }

          if (click_on_plus === true) {
            var old_scaleX = this.scaleX
            var old_scaleY = this.scaleY
            this.scaleX = this.scaleX*1.2;
            this.scaleY = this.scaleY*1.2;
            this.last_mouse1X = this.last_mouse1X - (width/(2*old_scaleX) - width/(2*this.scaleX));
            this.last_mouse1Y = this.last_mouse1Y - (height/(2*old_scaleY) - height/(2*this.scaleY));

          } else if (click_on_minus === true) {
            var old_scaleX = this.scaleX
            var old_scaleY = this.scaleY
            this.scaleX = this.scaleX/1.2;
            this.scaleY = this.scaleY/1.2;
            this.last_mouse1X = this.last_mouse1X - (width/(2*old_scaleX) - width/(2*this.scaleX));
            this.last_mouse1Y = this.last_mouse1Y - (height/(2*old_scaleY) - height/(2*this.scaleY));

          } else if (click_on_zoom_window === true) {
            this.zw_bool = !this.zw_bool;
            this.select_bool = false;

          } else if (click_on_reset === true){
            this.scaleX = this.init_scale;
            this.scaleY = this.init_scale;
            this.scale = this.init_scale;
            this.scroll_x = 0;
            this.scroll_y = 0;
            this.last_mouse1X = (width/2 - (coeff_pixel*this.maxX - coeff_pixel*this.minX)*this.scaleX/2)/this.scaleX - coeff_pixel*this.minX;
            this.last_mouse1Y = (height/2 - (coeff_pixel*this.maxY - coeff_pixel*this.minY)*this.scaleY/2)/this.scaleY - coeff_pixel*this.minY;

          } else if (click_on_select === true) {
            this.zw_bool = false;
            this.select_bool = !this.select_bool;

          } 

          this.draw(false, 0, this.last_mouse1X, this.last_mouse1Y, this.scaleX, this.scaleY, width, height, plot_datas);
          this.draw(true, 0, this.last_mouse1X, this.last_mouse1Y, this.scaleX, this.scaleY, width, height, plot_datas);
        }
      isDrawing = false;
      mouse_mouving = false;
    })

    canvas.addEventListener('wheel', e => {
      var zoom_coeff = 1.1;
      var event = -e.deltaY;
      mouse3X = e.offsetX;
      mouse3Y = e.offsetY;
      if ((mouse3Y>=height - 25) && (mouse3X>25)) {
          var old_scaleX = this.scaleX;
          if (event>0) {
            this.scaleX = this.scaleX*zoom_coeff;
          } else {
            this.scaleX = this.scaleX/zoom_coeff;
          }
          this.scroll_x = this.scroll_x - e.deltaY/Math.abs(e.deltaY);
          this.last_mouse1X = this.last_mouse1X - ((width/2)/old_scaleX - (width/2)/this.scaleX);
      } else if ((mouse3X<=25) && (mouse3Y<height - 25)) {
          var old_scaleY = this.scaleY;
          if (event>0) {
            this.scaleY = this.scaleY*zoom_coeff;
          } else {
            this.scaleY = this.scaleY/zoom_coeff;
          }
          this.scroll_y = this.scroll_y - e.deltaY/Math.abs(e.deltaY);
          this.last_mouse1Y = this.last_mouse1Y - ((height/2)/old_scaleY - (height/2)/this.scaleY);
      } else {
          var old_scaleY = this.scaleY;
          var old_scaleX = this.scaleX;
          if (event>0) {
            this.scaleX = this.scaleX*zoom_coeff;
            this.scaleY = this.scaleY*zoom_coeff;
          } else {
            this.scaleX = this.scaleX/zoom_coeff;
            this.scaleY = this.scaleY/zoom_coeff;
          }
          
          this.scroll_x = this.scroll_x - e.deltaY/Math.abs(e.deltaY);
          this.scroll_y = this.scroll_y - e.deltaY/Math.abs(e.deltaY);
          this.last_mouse1X = this.last_mouse1X - (mouse3X/old_scaleX - mouse3X/this.scaleX);
          this.last_mouse1Y = this.last_mouse1Y - (mouse3Y/old_scaleY - mouse3Y/this.scaleY);
      }
      this.draw(false, 0, this.last_mouse1X, this.last_mouse1Y, this.scaleX, this.scaleY, width, height, plot_datas);
      this.draw(true, 0, this.last_mouse1X, this.last_mouse1Y, this.scaleX, this.scaleY, width, height, plot_datas);
    })
  }

  contains_undefined(list) {
    for (var i=0; i<list.length; i++) {
      if (typeof list[i] === "undefined") {
        return true;
      }
    }
    return false;
  }

  remove_selection(val, list){
    var temp = [];
    for (var i = 0; i < list.length; i++) {
      var d = list[i];
      if (val != d) {
        temp.push(d);
      }
    }
    return temp;
  }

  is_include(val, list){
    for (var i = 0; i < list.length; i++) {
      var d = list[i];
      if (val == d) {
        return true;
      }
    }
    return false;
  }
}

export class PlotContour extends PlotData {
  constructor(public datas:any, 
                public width: number,
                public height: number,
                public coeff_pixel: number) {
    super();
    var plot_datas = [];
    for (var i = 0; i < datas.length; i++) {
      var d = plot_datas[i];
      var a = PlotDataContour2D.deserialize(d);
      this.minX = Math.min(this.minX, a.minX);
      this.maxX = Math.max(this.maxX, a.maxX);
      this.minY = Math.min(this.minY, a.minY);
      this.maxY = Math.max(this.maxY, a.maxY);
      this.colour_to_plot_data[a.mouse_selection_color] = a;
      plot_datas.push(a);
    }
    this.define_canvas(width, height);
    this.mouse_interaction(width, height, plot_datas, coeff_pixel);
  }


}

export class PlotScatter extends PlotData {
  plot_datas:any;

  constructor(public datas:any, 
    public width: number,
    public height: number,
    public coeff_pixel: number) {
      super();
      this.plot_datas = [];
      for (var i = 0; i < datas.length; i++) {
        var d = datas[i];
        if (d['type'] == 'point') {
          var a = PlotDataPoint2D.deserialize(d)
          this.minX = Math.min(this.minX, a.minX);
          this.maxX = Math.max(this.maxX, a.maxX);
          this.minY = Math.min(this.minY, a.minY);
          this.maxY = Math.max(this.maxY, a.maxY);
          this.colour_to_plot_data[a.mouse_selection_color] = a;
          this.plot_datas.push(a);
        } else if (d['type'] == 'plot') {
          var b = PlotDataScatterPlot.deserialize(d);
          this.plot_datas.push(b);
        }
      }
      this.define_canvas(width, height);
      this.mouse_interaction(width, height, this.plot_datas, coeff_pixel);
  }

  draw_initial() {
    super.draw_initial(this.width, this.height, this.coeff_pixel);
    this.draw(true, 0, this.last_mouse1X, this.last_mouse1Y, this.scaleX, this.scaleY, this.width, this.height, this.plot_datas);
    this.draw(false, 0, this.last_mouse1X, this.last_mouse1Y, this.scaleX, this.scaleY, this.width, this.height, this.plot_datas);
  }
  
}

class MyMath {
  public static round(x:number, n:number) {
    return Math.round(x*Math.pow(10,n)) / Math.pow(10,n);
  }
}

class Shape {
  public static crux(context:any, cx:number, cy:number, length:number) {
    context.moveTo(cx, cy);
    context.lineTo(cx - length, cy);
    context.moveTo(cx, cy);
    context.lineTo(cx + length, cy);
    context.moveTo(cx, cy);
    context.lineTo(cx, cy - length);
    context.moveTo(cx, cy);
    context.lineTo(cx, cy + length);
  }

  public static roundRect(x, y, w, h, radius, context) {
    var r = x + w;
    var b = y + h;
    context.beginPath();
    context.strokeStyle="black";
    context.lineWidth="3";
    context.moveTo(x+radius, y);
    context.lineTo(r-radius, y);
    context.quadraticCurveTo(r, y, r, y+radius);
    context.lineTo(r, y+h-radius);
    context.quadraticCurveTo(r, b, r-radius, b);
    context.lineTo(x+radius, b);
    context.quadraticCurveTo(x, b, x, b-radius);
    context.lineTo(x, y+radius);
    context.quadraticCurveTo(x, y, x+radius, y);
    context.stroke();
    context.closePath();
  }

  public static Is_in_rect(x, y, rect_x, rect_y, rect_w, rect_h) {
    return ((x>=rect_x) && (x<= rect_x + rect_w) && (y>=rect_y) && (y<=rect_y + rect_h))
  }

  public static createButton(x, y, w, h, context, text) {
    context.beginPath();
    context.fillStyle = 'white';
    context.lineWidth = "3";
    context.rect(x,y,w,h);
    context.stroke();
    context.fill();
    context.closePath();
    context.beginPath();
    context.fillStyle = "black"
    context.textAlign = "center";
    context.font = "12px Arial";
    context.fillText(text, x+w/2, y+h/1.8);
    context.fill();
    context.closePath();
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
        this.minX = Math.min(this.minX, d.minX);
        this.maxX = Math.max(this.maxX, d.maxX);
        this.minY = Math.min(this.minY, d.minY);
        this.maxY = Math.max(this.maxY, d.maxY);
      }
      this.mouse_selection_color = genColor();

  }

  public static deserialize(serialized) {
      var temp = serialized['plot_data_states'];
      var plot_data_states = [];
      for (var i = 0; i < temp.length; i++) {
        var d = temp[i];
        plot_data_states.push(PlotDataState.deserialize(d));
      }
      var temp = serialized['plot_data_primitives'];
      var plot_data_primitives = [];
      for (var i = 0; i < temp.length; i++) {
        var d = temp[i];
        if (d['type'] == 'line') {
          plot_data_primitives.push(PlotDataLine2D.deserialize(d));
        }
        if (d['type'] == 'circle') {
          plot_data_primitives.push(PlotDataCircle2D.deserialize(d));
        }
        if (d['type'] == 'arc') {
          plot_data_primitives.push(PlotDataArc2D.deserialize(d));
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
      var temp = serialized['plot_data_states'];
      var plot_data_states = [];
      for (var i = 0; i < temp.length; i++) {
        var d = temp[i];
        plot_data_states.push(PlotDataState.deserialize(d));
      }
      return new PlotDataLine2D(serialized['data'],
                               plot_data_states,
                               serialized['type'],
                               serialized['name']);
  }

  draw(context, first_elem, mvx, mvy, scaleX, scaleY) {
    if (first_elem) {
      context.moveTo(scaleX*(1000*this.data[0]+ mvx), scaleX*(1000*this.data[1]+ mvy));
    }
    context.lineTo(scaleY*(1000*this.data[2]+ mvx), scaleY*(1000*this.data[3]+ mvy));
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

  draw(context, first_elem, mvx, mvy, scaleX, scaleY) {
    context.arc(scaleX*(1000*this.cx+ mvx), scaleY*(1000*this.cy+ mvy), scaleX*1000*this.r, 0, 2*Math.PI);
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
      var plot = this.plot_data_states[i];
      var point_size = plot.point_size.size;
      if (point_size==1||point_size==2||point_size==3||point_size==4) {
        var height = plot.window_size.height;
        var width = plot.window_size.width;
      } else {
        throw new Error('Invalid point_size');
      }
    }
    this.size = point_size * Math.min(height,width)/150;
    this.minX = this.cx - this.size;
    this.maxX = this.cx + this.size;
    this.minY = this.cy - this.size;
    this.maxY = this.cy + this.size;
    this.mouse_selection_color = genColor();
    }

    public static deserialize(serialized) {
      var temp = serialized['plot_data_states'];
      var plot_data_states = [];
      for (var i = 0; i < temp.length; i++) {
        var d = temp[i];
        plot_data_states.push(PlotDataState.deserialize(d));
      }
      return new PlotDataPoint2D(serialized['data'],
                                  serialized['cx'],
                                  serialized['cy'],
                                  serialized['size'],
                                  plot_data_states,
                                  serialized['type'],
                                  serialized['name']);
    }

    draw(context, context_hidden, mvx, mvy, scaleX, scaleY, init_scale) {
        for (var i=0; i<this.plot_data_states.length; i++) {
          context.lineWidth = this.plot_data_states[i].stroke_width;
          var shape = this.plot_data_states[i].shape_set.shape;
          if (shape == 'circle') {
            context.arc(scaleX*(1000*this.cx+ mvx), scaleY*(1000*this.cy+ mvy), init_scale*1000*this.size, 0, 2*Math.PI);
            context.stroke();
          } else if (shape == 'square') {
            context.rect(scaleX*(1000*(this.cx - this.size) + mvx),scaleY*(1000*(this.cy - this.size) + mvy),init_scale*1000*this.size*2, init_scale*1000*this.size*2);
            context.stroke();
          } else if (shape == 'crux') {
            context.rect(scaleX*(1000*this.cx + mvx), scaleY*(1000*this.cy + mvy),init_scale*1000*this.size, init_scale*100*this.size);
            context.rect(scaleX*(1000*this.cx + mvx), scaleY*(1000*this.cy + mvy),-init_scale*1000*this.size, init_scale*100*this.size);
            context.rect(scaleX*(1000*this.cx + mvx), scaleY*(1000*this.cy + mvy),init_scale*100*this.size, init_scale*1000*this.size);
            context.rect(scaleX*(1000*this.cx + mvx), scaleY*(1000*this.cy + mvy),init_scale*100*this.size, -init_scale*1000*this.size);
            context.fillStyle = context.strokeStyle;
            context.stroke();

          } else {
            throw new Error('Invalid shape for point');
          }
        }
    }
}

export class PlotDataScatterPlot {
  colorStroke:any;
  x_step:number;
  y_step:number;
  constructor(public nb_points_x:number,
                     public nb_points_y:number,
                     public font_size:number,
                     public graduation_color:string,
                     public name:string, 
                     public type:string, 
                     public plot_data_states:PlotDataState[]) {

    for (var i=0; i<this.plot_data_states.length; i++) {
      var plot = this.plot_data_states[i];
      this.colorStroke = plot.color_line;
    }
  }

  public static deserialize(serialized) {
    var temp = serialized['plot_data_states'];
    var plot_data_states = [];
    for (var i = 0; i < temp.length; i++) {
      var d = temp[i];
      plot_data_states.push(PlotDataState.deserialize(d));
    }
    return new PlotDataScatterPlot(serialized['nb_points_x'],
                                  serialized['nb_points_y'],
                                  serialized['font_size'],
                                  serialized['graduation_color'],
                                  serialized['name'],
                                  serialized['type'],
                                  serialized['plot_data_states']);
  }

  draw_axes(context, mvx, mvy, scaleX, scaleY, width, height, minX, maxX, minY, maxY, x_step, y_step) {
    //pour l'axe des x
    var i=0;
    context.textAlign = 'center';
    var x_nb_digits = Math.max(0, 1-Math.floor(Math.log10(x_step)));
    while(minX + i*x_step < maxX) {
      if ((scaleX*(1000*(minX + i*x_step) + mvx) >0) && (scaleX*(1000*(minX + i*x_step) + mvx)<width)) {
        context.moveTo(scaleX*(1000*(minX + i*x_step) + mvx), height - 23);
        context.lineTo(scaleX*(1000*(minX + i*x_step) + mvx), height - 17);
        context.fillText(MyMath.round(minX + i*x_step, x_nb_digits), scaleX*(1000*(minX + i*x_step) + mvx), height - 4 );
      } 
      i++
    }
    if ((scaleX*(1000*(minX + i*x_step) + mvx) >0) && (scaleX*(1000*(minX + i*x_step) + mvx)<width)) {
      context.moveTo(scaleX*(1000*(minX + i*x_step) + mvx), height - 23);
      context.lineTo(scaleX*(1000*(minX + i*x_step) + mvx), height - 17);
      context.fillText(MyMath.round(minX + i*x_step, x_nb_digits), scaleX*(1000*(minX + i*x_step) + mvx), height - 4 );
    }
    
      //pour l'axe des y
    i=0
    var real_minY = -maxY;
    var real_maxY = -minY;
    context.textAlign = 'start';
    var y_nb_digits = Math.max(0, 1-Math.floor(Math.log10(y_step)));
    while (real_minY + (i-1)*y_step < real_maxY) {
      if ((scaleY*(-1000*(real_minY + i*y_step) + mvy)>0) && (scaleY*(-1000*(real_minY + i*y_step) + mvy)<height)) {
        context.moveTo(7, scaleY*(-1000*(real_minY + i*y_step) + mvy));
        context.lineTo(13, scaleY*(-1000*(real_minY + i*y_step) + mvy));
        context.fillText(MyMath.round(real_minY + i*y_step, y_nb_digits), 15, scaleY*(-1000*(real_minY + i*y_step) + mvy) + 5);
      }
      i++;
    }
    if ((scaleY*(-1000*(real_minY + i*y_step) + mvy)>0) && (scaleY*(-1000*(real_minY + i*y_step) + mvy)<height)) {
      context.moveTo(7, scaleY*(-1000*(real_minY + i*y_step) + mvy));
      context.lineTo(13, scaleY*(-1000*(real_minY + i*y_step) + mvy));
      context.fillText(MyMath.round(real_minY + i*y_step, y_nb_digits), 15, scaleY*(-1000*(real_minY + i*y_step) + mvy) + 5);
    }
    context.stroke();
  }

  draw(context, mvx, mvy, scaleX, scaleY, width, height, init_scale, minX, maxX, minY, maxY, scroll_x, scroll_y) {
    // Dessin du repère
    context.strokeStyle = this.colorStroke;

    //Flèches
    context.moveTo(0, 20);
    context.lineTo(10, 0);
    context.moveTo(10, 0);
    context.lineTo(20, 20);
    
    context.moveTo(width - 20, height - 30);
    context.lineTo(width, height - 20);
    context.moveTo(width, height - 20);
    context.lineTo(width - 20, height - 10);

    //Axes
    context.moveTo(10, 0);
    context.lineTo(10, height);

    context.moveTo(0, height - 20);
    context.lineTo(width, height - 20);
    //Graduations

    if (scroll_x % 5 == 0) {
      var kx = scaleX/init_scale;
      this.x_step = (maxX - minX)/(kx*(this.nb_points_x-1));
    } 
    if (scroll_y % 5 == 0) {
      var ky = scaleY/init_scale;
      this.y_step = (maxY - minY)/(ky*(this.nb_points_y-1));
    }

    context.font = this.font_size.toString() + 'px Arial';
    context.fillStyle = this.graduation_color;
    
    this.draw_axes(context, mvx, mvy, scaleX, scaleY, width, height, minX, maxX, minY, maxY, this.x_step, this.y_step);
    
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
      var plot_data_states = [];
      for (var i = 0; i < temp.length; i++) {
        var d = temp[i];
        plot_data_states.push(PlotDataState.deserialize(d));
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

  draw(context, first_elem, mvx, mvy, scaleX, scaleY) {
    var ptsa = [];
    for (var l = 0; l < this.data.length; l++) {
      ptsa.push(scaleX*(1000*this.data[l]['x']+ mvx));
      ptsa.push(scaleY*(1000*this.data[l]['y']+ mvy));
    }
    var tension = 0.4;
    var isClosed = false;
    var numOfSegments = 16;
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
        color_surface = ColorSurfaceSet.deserialize(serialized['color_surface']);
      }
      var hatching = null;
      if (serialized['hatching'] != null) {
        hatching = HatchingSet.deserialize(serialized['hatching']);
      }
      var shape_set = null;
      if (serialized['shape_set'] != null) {
        shape_set = PointShapeSet.deserialize(serialized['shape_set']);
      }
      var window_size = null;
      if(serialized['window_size'] != null) {
        window_size = WindowSizeSet.deserialize(serialized['window_size']);
      }
      var point_size = null;
      if (serialized['point_size'] != null) {
        point_size = PointSizeSet.deserialize(serialized['point_size']);
      }
      var point_color = null;
      if (serialized['point_color'] != null) {
        point_color = PointColorSet.deserialize(serialized['point_color']);
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
    var nb_hatch = 20;
    var max_size = nb_hatch*this.hatch_spacing;

    var p_hatch = document.createElement("canvas");
    p_hatch.width = max_size;
    p_hatch.height = max_size;
    var pctx = p_hatch.getContext("2d");
    pctx.lineCap = 'square';
    pctx.strokeStyle = 'black';
    pctx.lineWidth = this.stroke_width;
    pctx.beginPath();
    var pos_x = - Math.pow(Math.pow(max_size,2)/2, 0.5);
    var pos_y = Math.pow(Math.pow(max_size,2)/2, 0.5);
    for (var i = 0; i <= 2*nb_hatch; i++) {
      pos_x = pos_x + this.hatch_spacing;
      pos_y = pos_y - this.hatch_spacing;
      pctx.moveTo(pos_x, pos_y);
      pctx.lineTo(pos_x + max_size, pos_y + max_size);
    }
    pctx.stroke();
    return p_hatch;
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

    nextCol += 50;
  }
  var col = "rgb(" + ret.join(',') + ")";
  return col;
}