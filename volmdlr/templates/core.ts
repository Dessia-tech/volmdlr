export class PlotData {
  context_show:any;
  context_hidden:any;
  minX:number;
  maxX:number;
  minY:number;
  maxY:number;
  init_scale:number;
  init_scaleX:number;
  init_scaleY:number;
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
  color_surface_on_mouse:string='lightskyblue';
  color_surface_on_click:string='blue';
  context:any;
  tooltip_ON:boolean = false;
  axis_ON:boolean = false;
  link_object_ON:boolean = false;
  index_first_in:number;
  index_last_in:number;
  nb_points_in:number;
  graph_ON:boolean=false;


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
    this.init_scale = Math.min(this.width/(this.coeff_pixel*this.maxX - this.coeff_pixel*this.minX), this.height/(this.coeff_pixel*this.maxY - this.coeff_pixel*this.minY));
    this.scale = this.init_scale;
    if ((this.axis_ON) && !(this.graph_ON)) {
      this.init_scaleX = (this.width-this.decalage_axis_x)/(this.coeff_pixel*this.maxX - this.coeff_pixel*this.minX);
      this.init_scaleY = (this.height - this.decalage_axis_y)/(this.coeff_pixel*this.maxY - this.coeff_pixel*this.minY);
      this.scaleX = this.init_scaleX;
      this.scaleY = this.init_scaleY;
      this.last_mouse1X = (this.width/2 - (this.coeff_pixel*this.maxX - this.coeff_pixel*this.minX)*this.scaleX/2)/this.scaleX - this.coeff_pixel*this.minX + this.decalage_axis_x/(2*this.scaleX);
      this.last_mouse1Y = (this.height/2 - (this.coeff_pixel*this.maxY - this.coeff_pixel*this.minY)*this.scaleY/2)/this.scaleY - this.coeff_pixel*this.minY - this.decalage_axis_y/(2*this.scaleY);
    } else if ((this.axis_ON) && (this.graph_ON)) {
      this.init_scaleX = (this.width-this.decalage_axis_x)/(this.coeff_pixel*this.maxX - this.coeff_pixel*this.minX);
      this.init_scaleY = (this.height - this.decalage_axis_y - (this.graph1_button_y + this.graph1_button_h + 5))/(this.coeff_pixel*this.maxY - this.coeff_pixel*this.minY);
      this.scaleX = this.init_scaleX;
      this.scaleY = this.init_scaleY;
      this.last_mouse1X = (this.width/2 - (this.coeff_pixel*this.maxX - this.coeff_pixel*this.minX)*this.scaleX/2)/this.scaleX - this.coeff_pixel*this.minX + this.decalage_axis_x/(2*this.scaleX);
      this.last_mouse1Y = (this.height/2 - (this.coeff_pixel*this.maxY - this.coeff_pixel*this.minY)*this.scaleY/2)/this.scaleY - this.coeff_pixel*this.minY - (this.decalage_axis_y - (this.graph1_button_y + this.graph1_button_h + 5))/(2*this.scaleY);
    } else {
      this.scaleX = this.init_scale;
      this.scaleY = this.init_scale;
      this.last_mouse1X = (this.width/2 - (this.coeff_pixel*this.maxX - this.coeff_pixel*this.minX)*this.scaleX/2)/this.scaleX - this.coeff_pixel*this.minX;
      this.last_mouse1Y = (this.height/2 - (this.coeff_pixel*this.maxY - this.coeff_pixel*this.minY)*this.scaleY/2)/this.scaleY - this.coeff_pixel*this.minY;
    }
    this.draw(true, 0, this.last_mouse1X, this.last_mouse1Y, this.scaleX, this.scaleY);
    this.draw(false, 0, this.last_mouse1X, this.last_mouse1Y, this.scaleX, this.scaleY);

  }

  draw_empty_canvas(hidden) {
    if (hidden) {
      this.context = this.context_hidden;
    } else {
      this.context = this.context_show;
    }
    this.context.clearRect(0, 0, this.width, this.height);
  }

  draw_contour(hidden, show_state, mvx, mvy, scaleX, scaleY, d) {
    if (d['type'] == 'contour') {
      this.context.beginPath();
      if (hidden) {
        this.context.fillStyle = d.mouse_selection_color;
      } else {
        this.context.strokeStyle = d.plot_data_states[show_state].color_line;
        this.context.lineWidth = d.plot_data_states[show_state].stroke_width;
        this.context.fillStyle = 'white';
        if (d.plot_data_states[show_state].hatching != null) {
          this.context.fillStyle = this.context.createPattern(d.plot_data_states[show_state].hatching.canvas_hatching,'repeat');
        }
        if (d.plot_data_states[show_state].color_surface != null) {
          this.context.fillStyle = d.plot_data_states[show_state].color_surface.color;
        }
        if (this.select_on_mouse == d) {
          this.context.fillStyle = this.color_surface_on_mouse;
        }
        for (var j = 0; j < this.select_on_click.length; j++) {
          var z = this.select_on_click[j];
          if (z == d) {
            this.context.fillStyle = this.color_surface_on_click;
          }
        }
      }
      for (var j = 0; j < d.plot_data_primitives.length; j++) {
        var elem = d.plot_data_primitives[j];
        if (j == 0) {var first_elem = true} else {var first_elem = false}
        elem.draw(this.context, first_elem,  mvx, mvy, scaleX, scaleY);
      }
      this.context.fill();
      this.context.stroke();
      this.context.closePath();
    
    }
  }

  draw_point(hidden, show_state, mvx, mvy, scaleX, scaleY, d) {
    if (d['type'] == 'point') {
      if (hidden) {
        this.context.fillStyle = d.mouse_selection_color;
      } else {
        this.context.fillStyle = d.plot_data_states[show_state].point_color.color_fill;
        this.context.lineWidth = d.plot_data_states[show_state].stroke_width;
        this.context.strokeStyle = d.plot_data_states[show_state].point_color.color_stroke;
        var shape = d.plot_data_states[show_state].shape_set.shape;

        if (shape == 'crux') {
          this.context.strokeStyle = d.plot_data_states[show_state].point_color.color_fill;
        }
        if (this.select_on_mouse == d) {
          this.context.fillStyle = this.color_surface_on_mouse;
        }
        for (var j = 0; j < this.select_on_click.length; j++) {
          var z = this.select_on_click[j];
          if (z == d) {
            if (shape == 'crux') {
              this.context.strokeStyle = this.color_surface_on_click;
            } else {
              this.context.fillStyle = this.color_surface_on_click;         
            }
          }
        }
      }
      var x = scaleX*(1000*d.cx+ mvx);
      var y = scaleY*(1000*d.cy + mvy);
      var length = 1000*d.size;

      var is_inside_canvas = ((x + length>=0) && (x - length <= this.width) && (y + length >= 0) && (y - length <= this.height));

      if (is_inside_canvas === true) {
        this.context.beginPath();
        d.draw(this.context, this.context_hidden, mvx, mvy, scaleX, scaleY);
        this.context.fill();
        this.context.stroke();
        this.context.closePath();
      }
    }
  }

  draw_axis(mvx, mvy, scaleX, scaleY, d) {
    if (d['type'] == 'axis'){
      this.context.beginPath();
      d.draw(this.context, mvx, mvy, scaleX, scaleY, this.width, this.height, this.init_scaleX, this.init_scaleY, this.minX, this.maxX, this.minY, this.maxY, this.scroll_x, this.scroll_y, this.decalage_axis_x, this.decalage_axis_y);
      this.context.closePath();
      this.context.fill();
    }
  }
  draw_tooltip(d, mvx, mvy) {
    if (d['type'] == 'tooltip') {
      this.tooltip_ON = true;
      d.manage_tooltip(this.context, mvx, mvy, this.scaleX, this.scaleY, this.init_scale, this.width, this.height, this.tooltip_list)
    }
  }

  draw_graph2D(d, hidden, mvx, mvy) {
    if ((d['type'] == 'graph2D') && (this.graph_to_display[d.id] === true)) {
      this.context.beginPath();
      this.context.setLineDash(d.dashline);
      this.context.strokeStyle = d.graph_colorstroke;
      this.context.lineWidth = d.graph_linewidth;
      for (var i=0; i<d.segments.length; i++) {
        if (i==0) {
          d.segments[i].draw(this.context, true, mvx, mvy, this.scaleX, this.scaleY);
        } else {
          d.segments[i].draw(this.context, false, mvx, mvy, this.scaleX, this.scaleY);
        }
      }
      this.context.stroke();
      this.context.setLineDash([]);

      [this.index_first_in, this.nb_points_in, this.index_last_in] = this.get_nb_points_inside_canvas(d.point_list, mvx, mvy);

      var step = d.display_step;
      if (this.nb_points_in<=10) {
        step = 1;
      }
      for (var i=0; i<d.point_list.length; i=i+step) {
        var point = d.point_list[i];
        this.draw_point(hidden, 0, mvx, mvy, this.scaleX, this.scaleY, point);
      }
    } else if ((d['type'] == 'graph2D') && (this.graph_to_display[d.id] === false)) {
      this.delete_clicked_points(d.point_list);
      this.delete_tooltip(d.point_list);
    }
  }

  draw_scatterplot(d, hidden, mvx, mvy) {
    if (d['type'] == 'ScatterPlot') {
      var new_point_list = this.refresh_point_list1(d.point_list,mvx,mvy);
      for (var i=0; i<new_point_list.length; i++) {
        var point = new_point_list[i];
        this.draw_point(hidden, 0, mvx, mvy, this.scaleX, this.scaleY, point);
      }
      // var to_unselect_list = [];
      // for (var i=0; i<this.select_on_click.length; i++) {
      //   if (!this.is_include(this.select_on_click[i],new_point_list)) {
      //     to_unselect_list.push(this.select_on_click[i]);
      //   }
      // }
      // this.delete_clicked_points(to_unselect_list);

      // var to_delete_tooltip_list = [];
      // for (var i=0; i<this.tooltip_list.length; i++) {
      //   if (!this.is_include(this.tooltip_list[i],new_point_list)) {
      //     to_delete_tooltip_list.push(this.tooltip_list[i]);
      //   }
      // }
      // this.delete_tooltip(to_delete_tooltip_list);
    } 
  }

  zoom_button(x, y, w, h) {
    if ((x<0) || (x+h>this.width) || (y<0) || (y+2*h>this.height)) {
      throw new Error("Invalid x or y, the zoom button is out of the canvas");
    }
    this.context.strokeStyle = 'black';
    this.context.beginPath();
    this.context.lineWidth = "2";
    this.context.fillStyle = 'white';
    this.context.rect(x, y, w, h);
    this.context.rect(x, y+h, w, h);
    this.context.moveTo(x, y+h);
    this.context.lineTo(x+w, y+h);
    Shape.crux(this.context, x+w/2, y+h/2, h/3);
    this.context.moveTo(x + w/2 - h/3, y + 3*h/2);
    this.context.lineTo(x + w/2 + h/3, y + 3*h/2);
    this.context.fill();
    this.context.stroke();
    this.context.closePath();
  }

  zoom_window_button(x, y, w, h) {
    if ((x<0) || (x+h>this.width) || (y<0) || (y+h>this.height)) {
      throw new Error("Invalid x or y, the zoom window button is out of the canvas");
    }
    this.context.strokeStyle = 'black';
    if (this.zw_bool) {
      Shape.createButton(x, y, w, h, this.context, "Z ON", "12px Arial");
    } else {
      Shape.createButton(x, y, w, h, this.context, "Z OFF", "12px Arial");
    }
    
  }

  reset_button(x, y, w, h) {
    if ((x<0) || (x+h>this.width) || (y<0) || (y+h>this.height)) {
      throw new Error("Invalid x or y, the reset button is out of the canvas");
    }
    this.context.strokeStyle = 'black';
    Shape.createButton(x, y, w, h, this.context, "Reset", "12px Arial");
  }

  selection_button(x, y, w, h) {
    if ((x<0) || (x+h>this.width) || (y<0) || (y+h>this.height)) {
      throw new Error("Invalid x or y, the selection button is out of the canvas");
    }
    this.context.strokeStyle = 'black';
    if (this.select_bool) {
      Shape.createButton(x, y, w, h, this.context, "S ON", "12px Arial")
    } else {
      Shape.createButton(x, y, w, h, this.context, "S OFF", "12px Arial")
    }
  }

  graph_buttons(x, y, w, h, police, text_spacing) {
    for (var i=0; i<this.nb_graph; i++) {
      this.context.font = police;
      var text = 'Graph ' + (i+1).toString();
      if (this.graph_to_display[i] === true) {
        Shape.createGraphButton(x + i*(w+text_spacing), y, w, h, this.context, text, police, this.graph_colorlist[i], false);
      } else {
        Shape.createGraphButton(x + i*(w+text_spacing), y, w, h, this.context, text, police, this.graph_colorlist[i], true);
      }
    }
  }

  mouse_interaction() {
    var isDrawing = false;
    var mouse_moving = false;
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
        mouse_moving = true;
        mouse2X = e.offsetX;
        mouse2Y = e.offsetY;
        this.draw(false, 0, this.last_mouse1X + mouse2X/this.scaleX - mouse1X/this.scaleX, this.last_mouse1Y + mouse2Y/this.scaleY - mouse1Y/this.scaleY, this.scaleX, this.scaleY);
        this.draw(true, 0, this.last_mouse1X + mouse2X/this.scaleX - mouse1X/this.scaleX, this.last_mouse1Y + mouse2Y/this.scaleY - mouse1Y/this.scaleY, this.scaleX, this.scaleY);
        
      } else if ((isDrawing === true) && (this.zw_bool||this.select_bool)) {
        mouse_moving = true;
        mouse2X = e.offsetX;
        mouse2Y = e.offsetY;
        this.draw(false, 0, this.last_mouse1X, this.last_mouse1Y, this.scaleX, this.scaleY);
        this.draw(true, 0, this.last_mouse1X, this.last_mouse1Y, this.scaleX, this.scaleY);
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
        this.draw(false, 0, this.last_mouse1X, this.last_mouse1Y, this.scaleX, this.scaleY);
        this.draw(true, 0, this.last_mouse1X, this.last_mouse1Y, this.scaleX, this.scaleY);
      }
    })

    canvas.addEventListener('mouseup', e => {

      var scale_ceil = 400*this.init_scale;
      var scale_floor = this.init_scale/3;

      var click_on_plus = Shape.Is_in_rect(mouse1X, mouse1Y, this.zoom_rect_x, this.zoom_rect_y, this.zoom_rect_w, this.zoom_rect_h);
      var click_on_minus = Shape.Is_in_rect(mouse1X, mouse1Y, this.zoom_rect_x, this.zoom_rect_y + this.zoom_rect_h, this.zoom_rect_w, this.zoom_rect_h);
      var click_on_zoom_window = Shape.Is_in_rect(mouse1X, mouse1Y, this.zw_x, this.zw_y, this.zw_w, this.zw_h);
      var click_on_reset = Shape.Is_in_rect(mouse1X, mouse1Y, this.reset_rect_x, this.reset_rect_y, this.reset_rect_w, this.reset_rect_h);
      var is_rect_big_enough = (Math.abs(mouse2X - mouse1X)>40) && (Math.abs(mouse2Y - mouse1Y)>30);
      var click_on_select = Shape.Is_in_rect(mouse1X, mouse1Y, this.select_x, this.select_y, this.select_w, this.select_h);
      var click_on_graph = false;
      for (var i=0; i<this.nb_graph; i++) {
        var click_on_graph_i = Shape.Is_in_rect(mouse1X, mouse1Y, this.graph1_button_x + i*(this.graph1_button_w + this.graph_text_spacing), this.graph1_button_y, this.graph1_button_w, this.graph1_button_h);
        click_on_graph = click_on_graph || click_on_graph_i;
      }
      var click_on_button = click_on_plus || click_on_minus || click_on_zoom_window || click_on_reset || click_on_select || click_on_graph;

      if (mouse_moving) {
          if ((this.zw_bool && is_rect_big_enough)) {
            var zoom_coeff_x = this.width/Math.abs(mouse2X - mouse1X);
            var zoom_coeff_y = this.height/Math.abs(mouse2Y - mouse1Y);
            if ((this.scaleX*zoom_coeff_x < scale_ceil) && (this.scaleY*zoom_coeff_y < scale_ceil)) {
              this.last_mouse1X = this.last_mouse1X - Math.min(mouse1X, mouse2X)/this.scaleX
              this.last_mouse1Y = this.last_mouse1Y - Math.min(mouse1Y,mouse2Y)/this.scaleY
              this.scaleX = this.scaleX*zoom_coeff_x;
              this.scaleY = this.scaleY*zoom_coeff_y;
              this.draw(false, 0, this.last_mouse1X, this.last_mouse1Y, this.scaleX, this.scaleY);
              this.draw(true, 0, this.last_mouse1X, this.last_mouse1Y, this.scaleX, this.scaleY);
            }
            
          } else if (this.select_bool) {
            for (var i=0; i<this.plot_datas.length; i++) {
              var d = this.plot_datas[i];
              var in_rect = Shape.Is_in_rect(this.scaleX*(1000*d.cx + this.last_mouse1X),this.scaleY*(1000*d.cy + this.last_mouse1Y), Math.min(mouse1X, mouse2X), Math.min(mouse1Y, mouse2Y), Math.abs(mouse2X - mouse1X), Math.abs(mouse2Y - mouse1Y));
              if ((d['type']=="point") && (in_rect === true) && !(this.is_include(d, this.select_on_click))) {
                this.select_on_click.push(d);
              } else if (d['type'] == 'graph2D') {
                for (var j=0; j<d.point_list.length; j++) {
                  var x = this.scaleX*(1000*d.point_list[j].cx + this.last_mouse1X);
                  var y = this.scaleY*(1000*d.point_list[j].cy + this.last_mouse1Y);
                  in_rect = Shape.Is_in_rect(x, y, Math.min(mouse1X, mouse2X), Math.min(mouse1Y, mouse2Y), Math.abs(mouse2X - mouse1X), Math.abs(mouse2Y - mouse1Y));
                  if ((in_rect===true) && !(this.is_include(d.point_list[j], this.select_on_click))) {
                    this.select_on_click.push(d.point_list[j])
                  }
                }
              }
            }
            this.draw(false, 0, this.last_mouse1X, this.last_mouse1Y, this.scaleX, this.scaleY);
            this.draw(true, 0, this.last_mouse1X, this.last_mouse1Y, this.scaleX, this.scaleY);

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
          if (this.tooltip_ON) {
              if (this.is_include(click_plot_data, this.tooltip_list) && (!this.is_include(click_plot_data, this.select_on_click))) {
                this.tooltip_list = this.remove_selection(click_plot_data, this.tooltip_list);
              } else if (!this.is_include(click_plot_data, this.tooltip_list) && this.is_include(click_plot_data, this.select_on_click)){
                this.tooltip_list.push(click_plot_data);
              }
          }
          
          if (this.contains_undefined(this.select_on_click) && !click_on_button) {
            this.select_on_click = [];
            this.tooltip_list = [];
          }
 
          if ((click_on_plus === true) && (this.scaleX*1.2 < scale_ceil) && (this.scaleY*1.2 < scale_ceil)) {
            var old_scaleX = this.scaleX
            var old_scaleY = this.scaleY
            this.scaleX = this.scaleX*1.2;
            this.scaleY = this.scaleY*1.2;
            this.last_mouse1X = this.last_mouse1X - (this.width/(2*old_scaleX) - this.width/(2*this.scaleX));
            this.last_mouse1Y = this.last_mouse1Y - (this.height/(2*old_scaleY) - this.height/(2*this.scaleY));
            this.scroll_x = 0;
            this.scroll_y = 0;

          } else if ((click_on_minus === true) && (this.scaleX/1.2 > scale_floor) && (this.scaleY/1.2 > scale_floor)) {
            var old_scaleX = this.scaleX
            var old_scaleY = this.scaleY
            this.scaleX = this.scaleX/1.2;
            this.scaleY = this.scaleY/1.2;
            this.last_mouse1X = this.last_mouse1X - (this.width/(2*old_scaleX) - this.width/(2*this.scaleX));
            this.last_mouse1Y = this.last_mouse1Y - (this.height/(2*old_scaleY) - this.height/(2*this.scaleY));
            this.scroll_x = 0;
            this.scroll_y = 0;

          } else if (click_on_zoom_window === true) {
            this.zw_bool = !this.zw_bool;
            this.select_bool = false;
            
          } else if (click_on_reset === true) {
            this.scaleX = this.init_scaleX;
            this.scaleY = this.init_scaleY;
            this.scale = this.init_scale;
            this.scroll_x = 0;
            this.scroll_y = 0;
            if ((this.axis_ON === true) && (this.graph_ON === false)) {
              this.last_mouse1X = (this.width/2 - (this.coeff_pixel*this.maxX - this.coeff_pixel*this.minX)*this.scaleX/2)/this.scaleX - this.coeff_pixel*this.minX + this.decalage_axis_x/(2*this.scaleX);
              this.last_mouse1Y = (this.height/2 - (this.coeff_pixel*this.maxY - this.coeff_pixel*this.minY)*this.scaleY/2)/this.scaleY - this.coeff_pixel*this.minY - this.decalage_axis_y/(2*this.scaleY);
            } else if ((this.axis_ON === true) && (this.graph_ON === true)) {
              this.last_mouse1X = (this.width/2 - (this.coeff_pixel*this.maxX - this.coeff_pixel*this.minX)*this.scaleX/2)/this.scaleX - this.coeff_pixel*this.minX + this.decalage_axis_x/(2*this.scaleX);
              this.last_mouse1Y = (this.height/2 - (this.coeff_pixel*this.maxY - this.coeff_pixel*this.minY)*this.scaleY/2)/this.scaleY - this.coeff_pixel*this.minY - (this.decalage_axis_y - (this.graph1_button_y + this.graph1_button_h + 5))/(2*this.scaleY);
            } else {
              this.last_mouse1X = (this.width/2 - (this.coeff_pixel*this.maxX - this.coeff_pixel*this.minX)*this.scaleX/2)/this.scaleX - this.coeff_pixel*this.minX;
              this.last_mouse1Y = (this.height/2 - (this.coeff_pixel*this.maxY - this.coeff_pixel*this.minY)*this.scaleY/2)/this.scaleY - this.coeff_pixel*this.minY;
            }
            
          } else if (click_on_select === true) {
            this.zw_bool = false;
            this.select_bool = !this.select_bool;

          } else if (click_on_graph) {
            for (var i=0; i<this.nb_graph; i++) {
              var click_on_graph_i = Shape.Is_in_rect(mouse1X, mouse1Y, this.graph1_button_x + i*(this.graph1_button_w + this.graph_text_spacing), this.graph1_button_y, this.graph1_button_w, this.graph1_button_h);
              if (click_on_graph_i === true) {
                this.graph_to_display[i] = !this.graph_to_display[i];
              }
            }
          }

          this.draw(false, 0, this.last_mouse1X, this.last_mouse1Y, this.scaleX, this.scaleY);
        }
      isDrawing = false;
      mouse_moving = false;
    })

    canvas.addEventListener('wheel', e => {
      var scale_ceil = 400*this.init_scale;
      var scale_floor = this.init_scale/100;
      var zoom_coeff = 1.1;
      var event = -e.deltaY;
      mouse3X = e.offsetX;
      mouse3Y = e.offsetY;
      if ((mouse3Y>=this.height - this.decalage_axis_y) && (mouse3X>this.decalage_axis_x) && this.axis_ON) {
          var old_scaleX = this.scaleX;
          if ((event>0) && (this.scaleX*zoom_coeff<scale_ceil)) {
            this.scaleX = this.scaleX*zoom_coeff;
            this.scroll_x = this.scroll_x - e.deltaY/Math.abs(e.deltaY);
            this.last_mouse1X = this.last_mouse1X - ((this.width/2)/old_scaleX - (this.width/2)/this.scaleX);
          } else if ((event<0) && this.scaleX/zoom_coeff>scale_floor) {
            this.scaleX = this.scaleX/zoom_coeff;
            this.scroll_x = this.scroll_x - e.deltaY/Math.abs(e.deltaY);
            this.last_mouse1X = this.last_mouse1X - ((this.width/2)/old_scaleX - (this.width/2)/this.scaleX);
          }         

      } else if ((mouse3X<=this.decalage_axis_x) && (mouse3Y<this.height - this.decalage_axis_y) && this.axis_ON) {
          var old_scaleY = this.scaleY;
          if ((event>0) && (this.scaleY*zoom_coeff<scale_ceil)) {
            this.scaleY = this.scaleY*zoom_coeff;
            this.scroll_y = this.scroll_y - e.deltaY/Math.abs(e.deltaY);
            this.last_mouse1Y = this.last_mouse1Y - ((this.height/2)/old_scaleY - (this.height/2)/this.scaleY);
          } else if ((event<0) && this.scaleY/zoom_coeff>scale_floor) {
            this.scaleY = this.scaleY/zoom_coeff;
            this.scroll_y = this.scroll_y - e.deltaY/Math.abs(e.deltaY);
            this.last_mouse1Y = this.last_mouse1Y - ((this.height/2)/old_scaleY - (this.height/2)/this.scaleY);
          }
          
      } else {
          var old_scaleY = this.scaleY;
          var old_scaleX = this.scaleX;
          if ((event>0) && (this.scaleX*zoom_coeff<scale_ceil) && (this.scaleY*zoom_coeff<scale_ceil)) {
            this.scaleX = this.scaleX*zoom_coeff;
            this.scaleY = this.scaleY*zoom_coeff;
            this.scroll_x = this.scroll_x - e.deltaY/Math.abs(e.deltaY);
            this.scroll_y = this.scroll_y - e.deltaY/Math.abs(e.deltaY);
            this.last_mouse1X = this.last_mouse1X - (mouse3X/old_scaleX - mouse3X/this.scaleX);
            this.last_mouse1Y = this.last_mouse1Y - (mouse3Y/old_scaleY - mouse3Y/this.scaleY);
          } else if ((event<0) && (this.scaleX/zoom_coeff>scale_floor) && (this.scaleY/zoom_coeff>scale_floor)) {
            this.scaleX = this.scaleX/zoom_coeff;
            this.scaleY = this.scaleY/zoom_coeff;
            this.scroll_x = this.scroll_x - e.deltaY/Math.abs(e.deltaY);
            this.scroll_y = this.scroll_y - e.deltaY/Math.abs(e.deltaY);
            this.last_mouse1X = this.last_mouse1X - (mouse3X/old_scaleX - mouse3X/this.scaleX);
            this.last_mouse1Y = this.last_mouse1Y - (mouse3Y/old_scaleY - mouse3Y/this.scaleY);
          }
        }
        this.draw(false, 0, this.last_mouse1X, this.last_mouse1Y, this.scaleX, this.scaleY);
        this.draw(true, 0, this.last_mouse1X, this.last_mouse1Y, this.scaleX, this.scaleY); 
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

  get_nb_points_inside_canvas(list_points, mvx, mvy) { //Sous hypothèse que la liste est ordonnée par ordre croissant des x
    var bool = true;
    var k = 0;
    var index_first_in = -1;
    var nb_points_in = 0;
    var index_last_in = -1;

    while ((k<list_points.length) && bool) {
      var x = this.scaleX*(1000*list_points[k].cx + mvx);
      var y = this.scaleY*(1000*list_points[k].cy + mvy);
      var is_inside_canvas = (x>=0) && (x<=this.width) && (y>=0) && (y<=this.height);
      if (is_inside_canvas === true) {
        index_first_in = k;
        bool = false;
      } else {
        k++;
      }
    }
    if (index_first_in == -1) {
      return [index_first_in, nb_points_in, index_last_in];
    }

    while (k<list_points.length) {
      var x = this.scaleX*(1000*list_points[k].cx + mvx);
      var y = this.scaleY*(1000*list_points[k].cy + mvy);
      var is_inside_canvas = (x>=0) && (x<=this.width) && (y>=0) && (y<=this.height);
      if (is_inside_canvas === true) {
        index_last_in = k;
        nb_points_in++;
      } 
      k++;
    }
    return [index_first_in, nb_points_in, index_last_in];
  }

  is_inside_canvas(point, mvx, mvy) {
    var x = this.scaleX*(1000*point.cx + mvx);
    var y = this.scaleY*(1000*point.cy + mvy);
    length = 100*point.size;
    return (x+length>=0) && (x<=this.width-length) && (y+length>=0) && (y-length<=this.height);
  }

  get_points_inside_canvas(list_points, mvx, mvy) {
    var new_list_points = [];
    for (var i=0; i<list_points.length; i++) {
      if (this.is_inside_canvas(list_points[i],mvx,mvy)) {
        new_list_points.push(list_points[i]);
      }
    }
    return new_list_points;
  }

  distance(p1,p2) {
    return Math.sqrt(Math.pow(p1[0] - p2[0],2) + Math.pow(p1[1] - p2[1], 2));
  }

  copy_list(list) {
    var new_list = [];
    for (var i=0; i<list.length; i++) {
      new_list.push(list[i]);
    }
    return new_list;
  }

  hashing_point(point, nb_x, nb_y, mvx, mvy) {
    var x_step = this.width/nb_x;
    var y_step = this.height/nb_y;
    var x = this.scaleX*(1000*point.cx + mvx);
    var y = this.scaleY*(1000*point.cy + mvy);
    var i = Math.ceil(x/x_step);
    var j = Math.ceil(y/y_step);
    var key = 100*i + j;
    return [key,point];
  }

  hashing_list(point_list, nb_x, nb_y, mvx, mvy) {
    var dict = [];
    for (var k=0; k<point_list.length; k++) {
      var point_dict = this.hashing_point(point_list[k], nb_x, nb_y, mvx, mvy);
      dict.push(point_dict);
    }
    return dict;
  }

  dehashing_list(hashed_point_list) {
    var point_list = [];
    for (var i=0; i<hashed_point_list.length; i++) {
      point_list.push(hashed_point_list[i][1]);
    }
    return point_list;
  }

  refresh_point_list1(point_list, mvx, mvy) {
    var new_point_list = this.get_points_inside_canvas(this.copy_list(point_list), mvx, mvy);
    var nb_x = 10;
    var nb_y = 10;
    var dict_point_list = this.hashing_list(new_point_list, nb_x, nb_y, mvx, mvy);
    var i=0;
    var length = dict_point_list.length;
    while (i<length) {
      var point_dict_i = dict_point_list[i];
      var i_0 = Math.floor(point_dict_i[0]/100);
      var i_1 = point_dict_i[0]%100; 
      var size_i = point_dict_i[1].size;
      var xi = this.scaleX*(1000*point_dict_i[1].cx + mvx);
      var yi = this.scaleY*(1000*point_dict_i[1].cy + mvy);
      var bool = false;
      var j = i+1;
      while (j<length) {
        var point_dict_j = dict_point_list[j];
        var j_0 = Math.floor(point_dict_j[0]/100);
        var j_1 = point_dict_j[0]%100;
        var size_j = point_dict_j[1].size;
        if (size_i>=size_j) {
          var max_size_index = i
          var min_size_index = j;
        } else {
          var max_size_index = j;
          var min_size_index = i;
        }
        var xj = this.scaleX*(1000*point_dict_j[1].cx + mvx);
        var yj = this.scaleY*(1000*point_dict_j[1].cy + mvy);
        var is_touching_ij = this.distance([xi,yi], [xj,yj])<1000*(size_i + size_j);
        if ((Math.abs(i_0-j_0)<=1) && (Math.abs(i_1-j_1)<=1) && is_touching_ij) {
          var copy_point_max_index = dict_point_list[max_size_index][1].copy();
          var copy_point_min_index = dict_point_list[min_size_index][1].copy();
          var new_cx = (copy_point_max_index.cx + copy_point_min_index.cx)/2;
          var new_cy = (copy_point_max_index.cy + copy_point_min_index.cy)/2;
          var copy_plot_data_states = [copy_point_max_index.plot_data_states[0]];
          var point = new PlotDataPoint2D([],new_cx, new_cy, copy_plot_data_states, 'point', '');
          var size_coeff = 1.15;
          point.size = dict_point_list[max_size_index][1].size*size_coeff;
          var point_i = dict_point_list[i][1];
          var point_j = dict_point_list[j][1];
          this.delete_clicked_points([point_i, point_j]);
          this.delete_tooltip([point_i, point_j]);
          dict_point_list = this.remove_selection(dict_point_list[i], dict_point_list);
          dict_point_list = this.remove_selection(dict_point_list[j-1], dict_point_list);
          dict_point_list.push(this.hashing_point(point, nb_x, nb_y, mvx, mvy));
          this.colour_to_plot_data[point.mouse_selection_color] = point;
          bool = true;
          break;
        } else {
          j++;
        }
      }
      if (bool) {length--} else {i++}
    }
    return this.dehashing_list(dict_point_list);
  }

  refresh_point_list(point_list, mvx, mvy) {
    var new_point_list = this.get_points_inside_canvas(this.copy_list(point_list), mvx, mvy);
    var i = 0;
    var length = new_point_list.length;
    while (i<length) {
      var size_i = new_point_list[i].size;
      var xi = this.scaleX*(1000*new_point_list[i].cx + mvx);
      var yi = this.scaleY*(1000*new_point_list[i].cy + mvy);
      var bool = false;
      var j = i+1;
      while (j<length) {
        var size_j = new_point_list[j].size;
        if (size_i>=size_j) {var max_size_index = i} else {var max_size_index = j}
        var xj = this.scaleX*(1000*new_point_list[j].cx + mvx);
        var yj = this.scaleY*(1000*new_point_list[j].cy + mvy);
        if (this.distance([xi,yi], [xj,yj])<1000*(new_point_list[i].size + new_point_list[j].size)) {
          var new_cx = (new_point_list[i].cx + new_point_list[j].cx)/2;
          var new_cy = (new_point_list[i].cy + new_point_list[j].cy)/2;
          var copy_plot_data_states = [new_point_list[max_size_index].plot_data_states[0].copy()];
          var point = new PlotDataPoint2D([],new_cx, new_cy, copy_plot_data_states, 'point', '');
          var size_coeff = 1.3;
          point.size = new_point_list[max_size_index].size*size_coeff;
          var point_i = new_point_list[i];
          var point_j = new_point_list[j];
          this.delete_clicked_points([point_i, point_j]);
          this.delete_tooltip([point_i, point_j]);
          new_point_list = this.remove_selection(new_point_list[i], new_point_list);
          new_point_list = this.remove_selection(new_point_list[j-1], new_point_list);
          new_point_list.push(point);
          this.colour_to_plot_data[point.mouse_selection_color] = point;
          bool = true;
          break;
        } else {
          j++
        }
      }
      if (bool) {
        length--;
      } else {
        i++;
      }
    }
    return new_point_list;
  }

  delete_clicked_points(point_list) {
    var i = 0;
    while (i<this.select_on_click.length) {
      if (this.is_include(this.select_on_click[i], point_list)) {
        this.select_on_click = this.remove_selection(this.select_on_click[i], this.select_on_click);
      } else {
        i++;
      }
    }
  }

  delete_tooltip(point_list) {
    var i = 0;
    while (i<this.tooltip_list.length) {
      if (this.is_include(this.tooltip_list[i], point_list)) {
        this.tooltip_list = this.remove_selection(this.tooltip_list[i], this.tooltip_list);
      } else {
        i++;
      }
    }
  }

}

export class PlotContour extends PlotData {
  plot_datas:any;
  constructor(public data:any, 
                public width: number,
                public height: number,
                public coeff_pixel: number) {
    super();
    this.plot_datas = [];
    for (var i = 0; i < data.length; i++) {
      var d = this.data[i];
      var a = PlotDataContour2D.deserialize(d);
      if (isNaN(this.minX)) {this.minX = a.minX} else {this.minX = Math.min(this.minX, a.minX)};
      if (isNaN(this.maxX)) {this.maxX = a.maxX} else {this.maxX = Math.max(this.maxX, a.maxX)};
      if (isNaN(this.minY)) {this.minY = a.minY} else {this.minY = Math.min(this.minY, a.minY)};
      if (isNaN(this.maxY)) {this.maxY = a.maxY} else {this.maxY = Math.max(this.maxY, a.maxY)};
      this.colour_to_plot_data[a.mouse_selection_color] = a;
      this.plot_datas.push(a);
    }
    this.define_canvas();
    this.mouse_interaction();
  }
  
  draw(hidden, show_state, mvx, mvy, scaleX, scaleY) {
    this.draw_empty_canvas(hidden);

    for (var i = 0; i < this.plot_datas.length; i++) {
      var d = this.plot_datas[i];
      this.draw_contour(hidden, show_state, mvx, mvy, scaleX, scaleY, d);
    }
  }
}

export class PlotScatter extends PlotData {
  plot_datas:any;
  tooltip_list:any[]=[];
  zoom_rect_x:number;
  zoom_rect_y:number;
  zoom_rect_w:number;
  zoom_rect_h:number;
  zw_bool:boolean;
  zw_x:number;
  zw_y:number;
  zw_w:number;
  zw_h:number;
  reset_rect_x:number;
  reset_rect_y:number;
  reset_rect_w:number;
  reset_rect_h:number;
  select_bool:boolean;
  select_x:number;
  select_y:number;
  select_w:number;
  select_h:number;
  sort_list_points:any[]=[];
  graph_to_display:boolean[]=[];
  graph1_button_x:number;
  graph1_button_y:number;
  graph1_button_w:number;
  graph1_button_h:number;
  nb_graph:number = 0;
  graph_colorlist:string[]=[];
  graph_text_spacing:number;
  decalage_axis_x = 50;
  decalage_axis_y = 20;
  last_point_list:any[]=[];

  constructor(public data:any, 
    public width: number,
    public height: number,
    public coeff_pixel: number) {
      super();
      this.zoom_rect_x = this.width - 45;
      this.zoom_rect_y = 10;
      this.zoom_rect_w = 35;
      this.zoom_rect_h = 25;
      this.zw_x = this.width - 45;
      this.zw_y = 70;
      this.zw_w = 35;
      this.zw_h = 30;
      this.reset_rect_x = this.width - 45;
      this.reset_rect_y = 110;
      this.reset_rect_w = 35;
      this.reset_rect_h = 30;
      this.select_x = this.width - 45;
      this.select_y = 150;
      this.select_w = 35;
      this.select_h = 30;
      this.graph1_button_y = 10;
      this.graph1_button_w = 30;
      this.graph1_button_h = 15;
      this.plot_datas = [];
      var graphID = 0;
      this.graph_text_spacing = 50;
      for (var i = 0; i < data.length; i++) {
        var d = data[i]; 
        var a;
        if (d['type'] == 'point') {
          a = PlotDataPoint2D.deserialize(d)
          if (isNaN(this.minX)) {this.minX = a.minX} else {this.minX = Math.min(this.minX, a.minX)};
          if (isNaN(this.maxX)) {this.maxX = a.maxX} else {this.maxX = Math.max(this.maxX, a.maxX)};
          if (isNaN(this.minY)) {this.minY = a.minY} else {this.minY = Math.min(this.minY, a.minY)};
          if (isNaN(this.maxY)) {this.maxY = a.maxY} else {this.maxY = Math.max(this.maxY, a.maxY)};
          this.colour_to_plot_data[a.mouse_selection_color] = a;
          this.plot_datas.push(a);
        
        } else if (d['type'] == 'axis') {
          this.axis_ON = true;
          a = PlotDataAxis.deserialize(d);
          this.plot_datas.push(a);
        } else if (d['type'] == 'tooltip') {
          a = PlotDataTooltip.deserialize(d);
          this.plot_datas.push(a);
        } else if (d['type'] == 'graph2D') {
          this.graph_ON = true;
          a = PlotDataGraph2D.deserialize(d);
          a.id = graphID;
          graphID++;
          this.graph_colorlist.push(a.point_list[0].plot_data_states[0].point_color.color_fill);
          this.graph_to_display.push(true);
          for (var j=0; j<a.point_list.length; j++) {
            var point = a.point_list[j];
            if (isNaN(this.minX)) {this.minX = point.minX} else {this.minX = Math.min(this.minX, point.minX)};
            if (isNaN(this.maxX)) {this.maxX = point.maxX} else {this.maxX = Math.max(this.maxX, point.maxX)};
            if (isNaN(this.minY)) {this.minY = point.minY} else {this.minY = Math.min(this.minY, point.minY)};
            if (isNaN(this.maxY)) {this.maxY = point.maxY} else {this.maxY = Math.max(this.maxY, point.maxY)};
            this.colour_to_plot_data[point.mouse_selection_color] = point;
          }
          this.plot_datas.push(a);
        } else if (d['type'] == 'ScatterPlot') {
          a = PlotDataScatter.deserialize(d);
          for (var j=0; j<a.point_list.length; j++) {
            var point = a.point_list[j];
            if (isNaN(this.minX)) {this.minX = point.minX} else {this.minX = Math.min(this.minX, point.minX)};
            if (isNaN(this.maxX)) {this.maxX = point.maxX} else {this.maxX = Math.max(this.maxX, point.maxX)};
            if (isNaN(this.minY)) {this.minY = point.minY} else {this.minY = Math.min(this.minY, point.minY)};
            if (isNaN(this.maxY)) {this.maxY = point.maxY} else {this.maxY = Math.max(this.maxY, point.maxY)};
            this.colour_to_plot_data[point.mouse_selection_color] = point;
          }
          this.plot_datas.push(a);
        }
      }
      this.nb_graph = graphID;
      this.graph1_button_x = width/2 - this.nb_graph*(this.graph1_button_w + this.graph_text_spacing)/2;
      this.define_canvas();
      this.mouse_interaction();
  }

  draw(hidden, show_state, mvx, mvy, scaleX, scaleY) {
    this.draw_empty_canvas(hidden);
    for (var i = 0; i < this.plot_datas.length; i++) {
      var d = this.plot_datas[i];
      this.draw_graph2D(d, hidden, mvx, mvy);
      this.draw_scatterplot(d, hidden, mvx, mvy);
      this.draw_point(hidden, show_state, mvx, mvy, scaleX, scaleY, d);
      this.draw_axis(mvx, mvy, scaleX, scaleY, d);
      this.draw_tooltip(d, mvx, mvy);
    }
      //Drawing the zooming button 
      this.zoom_button(this.zoom_rect_x, this.zoom_rect_y, this.zoom_rect_w, this.zoom_rect_h);
      
      //Drawing the button for zooming window selection
      this.zoom_window_button(this.zw_x,this.zw_y,this.zw_w,this.zw_h);
  
      //Drawing the reset button
      this.reset_button(this.reset_rect_x, this.reset_rect_y, this.reset_rect_w, this.reset_rect_h);
      
      //Drawing the selection button
      this.selection_button(this.select_x, this.select_y, this.select_w, this.select_h);

      //Drawing the enable/disable graph button
      this.graph_buttons(this.graph1_button_x ,this.graph1_button_y, this.graph1_button_w, this.graph1_button_h, '10px Arial', this.graph_text_spacing);
    
  }
}

class MyMath {
  public static round(x:number, n:number) {
    return Math.round(x*Math.pow(10,n)) / Math.pow(10,n);
  }
}

class Shape {

  public static drawLine(context, start, end) {
    context.moveTo(start[0], start[1]);
    context.lineTo(end[0], end[1]);
  }

  public static crux(context:any, cx:number, cy:number, length:number) {
    this.drawLine(context, [cx, cy], [cx - length, cy]);
    this.drawLine(context, [cx, cy], [cx + length, cy]);
    this.drawLine(context, [cx, cy], [cx, cy - length]);
    this.drawLine(context, [cx, cy], [cx, cy + length]);
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

  public static createButton(x, y, w, h, context, text, police) {
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
    context.font = police;
    context.fillText(text, x+w/2, y+h/1.8);
    context.fill();
    context.closePath();
  }

  public static createGraphButton(x, y, w, h, context, text, police, colorfill, strikeout) {
    context.beginPath();
    context.fillStyle = colorfill;
    context.rect(x,y,w,h);
    context.fill();
    context.closePath();
    context.beginPath();
    context.fillStyle = 'grey';
    context.textAlign = 'start';
    context.textBaseline = 'middle';
    context.font = police;
    context.fillText(text, x+w + 5, y+h/1.8);
    context.fill();
    if (strikeout === true) {
      var text_w = context.measureText(text).width;
      context.lineWidth = 1.5;
      context.strokeStyle = 'grey';
      Shape.drawLine(context, [x+w + 5, y+h/1.8], [x+w+5+text_w, y+h/2]);
      context.stroke();
    }
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
      context.moveTo(scaleX*(1000*this.data[0]+ mvx), scaleY*(1000*this.data[1]+ mvy));
    }
    context.lineTo(scaleX*(1000*this.data[2]+ mvx), scaleY*(1000*this.data[3]+ mvy));
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
  size:number;
  k:number=1;

  constructor(public data:any,
              public cx:number,
              public cy:number,
              public plot_data_states:PlotDataState[],
              public type:string,
              public name:string) {
    
    var plot = plot_data_states[0];
    var point_size = plot.point_size.size;
    if (point_size<1) {
      throw new Error('Invalid point_size');
    }
    this.size = this.k*point_size/400;
    this.minX = this.cx - 2.5*this.size;
    this.maxX = this.cx + 2.5*this.size;
    this.minY = this.cy - 5*this.size;
    this.maxY = this.cy + 5*this.size;
    
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
                                  plot_data_states,
                                  serialized['type'],
                                  serialized['name']);
    }

    draw(context, context_hidden, mvx, mvy, scaleX, scaleY) {
      var show_states = 0;
        var shape = this.plot_data_states[show_states].shape_set.shape;
        if (shape == 'circle') {
          context.arc(scaleX*(1000*this.cx+ mvx), scaleY*(1000*this.cy+ mvy), 1000*this.size, 0, 2*Math.PI);
          context.stroke();
        } else if (shape == 'square') {
          context.rect(scaleX*(1000*this.cx + mvx) - 1000*this.size,scaleY*(1000*this.cy + mvy) - 1000*this.size,1000*this.size*2, 1000*this.size*2);
          context.stroke();
        } else if (shape == 'crux') {
          context.rect(scaleX*(1000*this.cx + mvx), scaleY*(1000*this.cy + mvy),1000*this.size, 100*this.size);
          context.rect(scaleX*(1000*this.cx + mvx), scaleY*(1000*this.cy + mvy),-1000*this.size, 100*this.size);
          context.rect(scaleX*(1000*this.cx + mvx), scaleY*(1000*this.cy + mvy),100*this.size, 1000*this.size);
          context.rect(scaleX*(1000*this.cx + mvx), scaleY*(1000*this.cy + mvy),100*this.size, -1000*this.size);
          context.fillStyle = context.strokeStyle;
          context.stroke();
        } else {
          throw new Error('Invalid shape for point');
        }
        
    }

    copy() {
      return new PlotDataPoint2D(this.data, this.cx, this.cy, this.plot_data_states, this.type, this.name);
    }
}

export class PlotDataAxis {
  colorStroke:any;
  x_step:number;
  y_step:number;
  constructor(public nb_points_x:number,
                     public nb_points_y:number,
                     public font_size:number,
                     public graduation_color:string,
                     public axis_color:string,
                     public name:string, 
                     public arrow_on:boolean,
                     public axis_width:string,
                     public grid_on:boolean,
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
    return new PlotDataAxis(serialized['nb_points_x'],
                                  serialized['nb_points_y'],
                                  serialized['font_size'],
                                  serialized['graduation_color'],
                                  serialized['axis_color'],
                                  serialized['name'],
                                  serialized['arrow_on'],
                                  serialized['axis_width'],
                                  serialized['grid_on'], 
                                  serialized['type'],
                                  plot_data_states);
  }

  draw_graduations(context, mvx, mvy, scaleX, scaleY, axis_x_start, axis_x_end, axis_y_start, axis_y_end, minX, maxX, minY, maxY, x_step, y_step, font_size) {
    //pour l'axe des x
    var i=0;
    context.textAlign = 'center';
    var x_nb_digits = Math.max(0, 1-Math.floor(Math.log10(x_step)));
    var delta_x = maxX - minX;
    var grad_beg_x = minX - 10*delta_x;
    var grad_end_x = maxX + 10*delta_x;
    while(grad_beg_x + i*x_step < grad_end_x) {
      if ((scaleX*(1000*(grad_beg_x + i*x_step) + mvx) >axis_x_start) && (scaleX*(1000*(grad_beg_x + i*x_step) + mvx)<axis_x_end)) {
        
        if (this.grid_on === true) {
          context.strokeStyle = 'lightgrey';
          Shape.drawLine(context, [scaleX*(1000*(grad_beg_x + i*x_step) + mvx), axis_y_start], [scaleX*(1000*(grad_beg_x + i*x_step) + mvx), axis_y_end + 3]);
        } else {
          Shape.drawLine(context, [scaleX*(1000*(grad_beg_x + i*x_step) + mvx), axis_y_end - 3], [scaleX*(1000*(grad_beg_x + i*x_step) + mvx), axis_y_end + 3]);
        }
        context.fillText(MyMath.round(grad_beg_x + i*x_step, x_nb_digits), scaleX*(1000*(grad_beg_x + i*x_step) + mvx), axis_y_end + font_size );
      } 
      i++
    }
    
      //pour l'axe des y
    i=0
    var real_minY = -maxY;
    var real_maxY = -minY;
    var delta_y = maxY - minY;
    var grad_beg_y = real_minY - 10*delta_y;
    var grad_end_y = real_maxY + 10*delta_y;
    context.textAlign = 'end';
    context.textBaseline = 'middle';
    var y_nb_digits = Math.max(0, 1-Math.floor(Math.log10(y_step)));
    while (grad_beg_y + (i-1)*y_step < grad_end_y) {
      if ((scaleY*(-1000*(grad_beg_y + i*y_step) + mvy) > axis_y_start) && (scaleY*(-1000*(grad_beg_y + i*y_step) + mvy) < axis_y_end)) {
        if (this.grid_on === true) {
          context.strokeStyle = 'lightgrey';
          Shape.drawLine(context,[axis_x_start - 3, scaleY*(-1000*(grad_beg_y + i*y_step) + mvy)], [axis_x_end, scaleY*(-1000*(grad_beg_y + i*y_step) + mvy)]);
        } else {
          Shape.drawLine(context, [axis_x_start - 3, scaleY*(-1000*(grad_beg_y + i*y_step) + mvy)], [axis_x_start + 3, scaleY*(-1000*(grad_beg_y + i*y_step) + mvy)]);
        }   
        context.fillText(MyMath.round(grad_beg_y + i*y_step, y_nb_digits), axis_x_start - 5, scaleY*(-1000*(grad_beg_y + i*y_step) + mvy));
      }
      i++;
    }

    context.stroke();
  }

  draw(context, mvx, mvy, scaleX, scaleY, width, height, init_scaleX, init_scaleY, minX, maxX, minY, maxY, scroll_x, scroll_y, decalage_axis_x, decalage_axis_y) {
    // Dessin du repère
    context.beginPath();
    context.strokeStyle = this.axis_color;
    context.lineWidth = this.axis_width;
    var axis_x_start = decalage_axis_x;
    var axis_x_end = width;
    var axis_y_start = 0;
    var axis_y_end = height - decalage_axis_y;
    //Flèches
    if (this.arrow_on === true) {
      Shape.drawLine(context, [axis_x_start - 10, axis_y_start + 20], [axis_x_start, axis_y_start]);
      Shape.drawLine(context, [axis_x_start, axis_y_start], [axis_x_start + 10, axis_y_start + 20]);
      
      Shape.drawLine(context, [axis_x_end - 20, axis_y_end - 10], [axis_x_end, axis_y_end]);
      Shape.drawLine(context, [axis_x_end, axis_y_end], [axis_x_end - 20, axis_y_end + 10]);
    }
    
    //Axes
    Shape.drawLine(context, [axis_x_start, axis_y_start], [axis_x_start, axis_y_end]);
    Shape.drawLine(context, [axis_x_start, axis_y_end], [axis_x_end, axis_y_end]);

    context.stroke();

    //Graduations
    if (scroll_x % 5 == 0) {
      var kx = 1.1*scaleX/init_scaleX;
      this.x_step = (maxX - minX)/(kx*(this.nb_points_x-1));
    } 
    if (scroll_y % 5 == 0) {
      var ky = 1.1*scaleY/init_scaleY;
      this.y_step = (maxY - minY)/(ky*(this.nb_points_y-1));
    }
    
    context.font = this.font_size.toString() + 'px Arial';
    context.fillStyle = this.graduation_color;
    context.strokeStyle = this.axis_color;
    
    this.draw_graduations(context, mvx, mvy, scaleX, scaleY, axis_x_start, axis_x_end, axis_y_start, axis_y_end, minX, maxX, minY, maxY, this.x_step, this.y_step, this.font_size);
    context.closePath();
    
  }
}

export class PlotDataTooltip {
  constructor(public colorfill:string, public font:string, public tp_width:number, public tp_radius:any, public to_plot_list:any, public plot_data_states:PlotDataState[],public type:string, public name:string) {}

  public static deserialize(serialized) {
    var temp = serialized['plot_data_states']
      var plot_data_states = [];
      for (var i = 0; i < temp.length; i++) {
        var d = temp[i];
        plot_data_states.push(PlotDataState.deserialize(d));
      }
      return new PlotDataTooltip(serialized['colorfill'],
                                  serialized['font'],
                                  serialized['tp_width'],
                                  serialized['tp_radius'],
                                  serialized['to_plot_list'],
                                  plot_data_states,
                                  serialized['type'],
                                  serialized['name']);
  }

  draw(context, object, mvx, mvy, scaleX, scaleY, init_scale, canvas_width, canvas_height) {
    context.beginPath();
    var textfills = [];
    for (var i=0; i<this.to_plot_list.length; i++) {
      if (this.to_plot_list[i] == 'cx') {
        textfills.push('x : ' + MyMath.round(object.cx,4).toString());
      } else if (this.to_plot_list[i] == 'cy') {
        textfills.push('y : ' + MyMath.round(-object.cy, 4).toString());
      } else if (this.to_plot_list[i] == 'shape') {
        textfills.push('shape : ' + object.plot_data_states[0]['shape_set']['shape']);
      }
    }

    var font_size = Number(this.font.split('px')[0]);
    var tp_height = (textfills.length + 0.25)*font_size ;
    var cx = object.cx;
    var cy = object.cy;
    var point_size = object.plot_data_states[0].point_size.size;
    var decalage = 2.5*point_size + 5
    var tp_x = scaleX*(1000*cx + mvx) + decalage;
    var tp_y = scaleY*(1000*cy + mvy) - 1/2*tp_height;

    if (tp_x + this.tp_width  > canvas_width) {
      tp_x = scaleX*(1000*cx + mvx) - decalage - this.tp_width;
    }
    if (tp_y < 0) {
      tp_y = scaleY*(1000*cy + mvy);
    }
    if (tp_y + tp_height > canvas_height) {
      tp_y = scaleY*(1000*cy + mvy) - tp_height;
    }

    Shape.roundRect(tp_x, tp_y, this.tp_width, tp_height, this.tp_radius, context);
    context.strokeStyle = 'black';
    context.fillStyle = this.colorfill;
    context.stroke();
    context.fill();
    context.fillStyle = 'black';
    context.textAlign = 'center';
    context.textBaseline = 'Alphabetic';

    var x_middle = tp_x + 1/2*this.tp_width;
    context.font = this.font;

    var current_y = tp_y + 0.75*font_size;
    for (var i=0; i<textfills.length; i++) {
      context.fillText(textfills[i], x_middle, current_y);
      current_y = current_y + font_size;
    }
    context.closePath();
  }

  manage_tooltip(context, mvx, mvy, scaleX, scaleY, init_scale, canvas_width, canvas_height, tooltip_list) {
    for (var i=0; i<tooltip_list.length; i++) {
      if (!(typeof tooltip_list[i] === "undefined")) {
        this.draw(context, tooltip_list[i], mvx, mvy, scaleX, scaleY, init_scale, canvas_width, canvas_height);
      }
    }
  }
}

export class PlotDataGraph2D {
  id:number=0;
  constructor(public point_list:PlotDataPoint2D[],
              public dashline: number[],
              public graph_colorstroke: string,
              public graph_linewidth: number,
              public segments:PlotDataLine2D[],
              public display_step:number,
              public plot_data_states: PlotDataState[],
              public type: string,
              public name:string) {}
  
  public static deserialize(serialized) {
    var temp = serialized['plot_data_states'];
    var plot_data_states = [];
    for (var i = 0; i < temp.length; i++) {
      var d = temp[i];
      plot_data_states.push(PlotDataState.deserialize(d));
    }
    var point_list = [];
    temp = serialized['serialized_point_list'];
    for (var i=0; i<temp.length; i++) {
      var d = temp[i];
      point_list.push(PlotDataPoint2D.deserialize(d));
    }

    var segments = [];
    temp = serialized['serialized_segments'];
    for (i=0; i<temp.length; i++) {
      var d = temp[i];
      segments.push(PlotDataLine2D.deserialize(d));
    }
    return new PlotDataGraph2D(point_list,
                           serialized['dashline'],
                           serialized['graph_colorstroke'],
                           serialized['graph_linewidth'],
                           segments,
                           serialized['display_step'],
                           plot_data_states,
                           serialized['type'],
                           serialized['name']);
  }
}

export class PlotDataScatter {
  constructor(public point_list:PlotDataPoint2D[],
              public plot_data_states:PlotDataState[],
              public type:string,
              public name:string) {}
  
  public static deserialize(serialized) {
    var temp = serialized['plot_data_states'];
    var plot_data_states = [];
    for (var i = 0; i < temp.length; i++) {
      var d = temp[i];
      plot_data_states.push(PlotDataState.deserialize(d));
    }
    var point_list = [];
    temp = serialized['serialized_point_list'];
    for (var i=0; i<temp.length; i++) {
      var d = temp[i];
      point_list.push(PlotDataPoint2D.deserialize(d));
    }
    return new PlotDataScatter(point_list,
                               plot_data_states,
                               serialized['type'],
                               serialized['name']);
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
  copy() {
    return new PlotDataState(this.color_surface, this.color_map, this.hatching, this.opacity, this.dash, this.marker, this.color_line, this.shape_set, this.point_size, this.point_color, this.window_size, this.stroke_width, this.name);
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