export class D3Unpacker {

  constructor(D3Data, data_container){
    var scale_factor = 1; // Remove statement (Binding ? Arguments ?)
    var font_size = 1.5; // Remove statement (Binding ? Arguments ?)

    // Style definition of the SVG
    var svg2 = data_container.append("defs")
        .append("pattern")
        .attr("id","diagonal-stripe-1")
        .attr("width","1")
        .attr("height","1")
        .attr("patternUnits","userSpaceOnUse")
        .attr("patternTransform", "rotate(45)")
        .append("rect")
        .attr("width","0.1")
        .attr("height","100")
        .attr("transform","translate(0,0)")
        .attr("fill","black");
    var color_triangle_load = d3.rgb(0, 0, 256)
    data_container.selectAll("defs")
        .append("marker")
        .attr("id", "triangle_load")
        .attr("refX", 6)
        .attr("refY", 6)
        .attr('viewBox', "0 0 12 12")
        .attr("markerUnits","userSpaceOnUse")
        .attr("orient", "auto")
        .attr("fill", "blue")
        .append("path")
        .attr("d", "M 0 0 12 6 0 12 3 6")
        .attr('opacity', 0.7)

    data_container.selectAll("defs")
        .append("marker")
        .attr("id", "triangle_quote_end")
        .attr("refX", 12)
        .attr("refY", 6)
        .attr("markerWidth", 1)
        .attr("markerHeight", 2)
        .attr('viewBox', "0 0 12 12")
        .attr("markerUnits","userSpaceOnUse")
        .attr("orient", "auto")
        .append("path")
        .attr("d", "M 0 0 L 12 6 L 0 12")
        .style("fill", "none")
        .attr("stroke", "black")
    data_container.selectAll("defs")
        .append("marker")
        .attr("id", "triangle_quote_start")
        .attr("refX", 0)
        .attr("refY", 6)
        .attr("markerWidth", 1)
        .attr("markerHeight", 2)
        .attr('viewBox', "0 0 12 12")
        .attr("markerUnits","userSpaceOnUse")
        .attr("orient", "auto")
        .append("path")
        .attr("d", "M 12 0 L 0 6 L 12 12")
        .style("fill", "none")
        .attr("stroke", "black")

    D3Data.forEach(function(dict_component){
      // Trace quote
      if (dict_component.type == 'line') {
        var lines_data = [];
        lines_data.push(dict_component);
        lines_data.forEach(function(d){
          var stroke_color = d3.rgb(0*255, 1*255, 1*255);
          data_container.append("line")
            .attr('id', 'line_geom')
            .attr("x1", 1000*d.data[0])
            .attr("y1", 1000*d.data[1])
            .attr("x2", 1000*d.data[2])
            .attr("y2", 1000*d.data[3])
            .attr("stroke-width", d.size)
            .attr("stroke", d.color)
            .attr("stroke-width-init", d.size)
        })
      }
      else if (dict_component.type == "wire"){
        var paths_data = [];
        paths_data.push(dict_component);
        paths_data.forEach(function(d){
          var draw_data = []
          d.data.forEach(function(dd){
            draw_data.push([1000*dd['x'],1000*dd['y']])
          })
          var lineGenerator = d3.line();
          var pathString = lineGenerator(draw_data);
          data_container.append("path")
              .attr('d', pathString)
              .attr('id', 'path_geom')
              .attr("stroke-width", d.size)
              .attr("stroke", d.color)
              .attr("fill", "none")
              .attr("stroke-width-init", d.size)
        })
      }
      else if (dict_component.type == "circle"){
        var circles_data = [];
        circles_data.push(dict_component);
        circles_data.forEach(function(d){
          data_container.append("circle")
              .attr('r', 1000*d.r)
              .attr('cx', 1000*d.cx)
              .attr('cy', 1000*d.cy)
              .attr('id', 'circle_geom')
              // .attr("stroke-width", d.size)
              // .attr("stroke", 'black'')
              .attr("fill", d.color)
              .attr("pointer-events", 'all')
              //.on("mouseover", highlight)
              //.on("mouseout", deemphasize);
        })
      }
      else if (dict_component.type == "arc"){
        var arcs_data = [];
        arcs_data.push(dict_component);

        arcs_data.forEach(function(d){
          var draw_data = []
          d.data.forEach(function(dd){
            draw_data.push([1000*dd['x'],1000*dd['y']])
          })
          var lineGenerator = d3.line();
          lineGenerator.curve(d3.curveCardinal.tension(0));
          var pathString = lineGenerator(draw_data);
          data_container.append("path")
              .attr('d', pathString)
              .attr('id', 'path_geom')
              .attr("stroke-width", d.size)
              .attr("stroke", d.color)
              .attr("fill", "none")
              .attr("stroke-width-init", d.size)
        })
      }
      else if (dict_component.type == "text"){
        var text_data = [];
        text_data.push(dict_component);
        text_data.forEach(function(d){
          data_container.append("text")
            .attr("id", 'text_geom')
            .attr("transform", "translate(" + 1000*(d.x_label) + ","
                + 1000*d.y_label + ") rotate("+d.rot_label+")")
            .attr("font-size", d.font_size+'px')
            .attr("font-size-init", d.font_size)
            .attr('text-anchor', "middle")
            .attr('baseline-shift', d.baseline_shift+"ex")
            .attr("font-weight", 0.1)
            .attr("font-family", 'sans-serif')
            .attr("fill", 'black')
            .text(d.label);
        })
      }
      if (dict_component.type == 'quote') {
        var explore = dict_component['plot_data'];
        explore.forEach(function(d){
            if (d.marker == null){
              data_container.append("line")
                .attr('id', 'line_quote')
                .attr("x1", 1000*d.data[0])
                .attr("y1", 1000*d.data[1])
                .attr("x2", 1000*d.data[2])
                .attr("y2", 1000*d.data[3])
                .attr("stroke-width", d.size)
                .attr("stroke", 'black')
                .attr("fill", "none")
                .attr("stroke-linecap", "round")
                .attr("stroke-width-init", d.size)
              }
            else {
              data_container.append("line")
                .attr('id', 'line_quote')
                .attr("x1", 1000*d.data[0])
                .attr("y1", 1000*d.data[1])
                .attr("x2", 1000*d.data[2])
                .attr("y2", 1000*d.data[3])
                .attr("stroke-width", d.size)
                .attr("stroke", 'black')
                .attr("marker-end", "url(#triangle_quote_end)")
                .attr("marker-start", "url(#triangle_quote_start)")
                .attr("stroke-width-init", d.size)
              }
        })
        data_container.append("text")
            .attr("id", 'texte_quote')
            .attr("transform", "translate(" + 1000*(dict_component.x_label) + ","
                + 1000*dict_component.y_label + ") rotate("+dict_component.rot_label+")")
            .attr("dy", -0.3+"em")
            .attr("font-size", font_size)
            .attr('text-anchor', "middle")
            .text(dict_component.label);
      }

      // // Trace line
      // else if (dict_component.type == 'line') {
      //   var explore = dict_component['plot_data'];
      //   explore.forEach(function(d){
      //       data_container.append("line")
      //         .attr('id', 'line_constructor')
      //         .attr("x1", 1000*d.data[0])
      //         .attr("y1", 1000*d.data[1])
      //         .attr("x2", 1000*d.data[2])
      //         .attr("y2", 1000*d.data[3])
      //         .attr("stroke-width", d.stroke_width)
      //         .attr("stroke", 'black')
      //         .attr("fill", "none")
      //         .attr("stroke-linecap", "round");
      //   })
      // }

      // Trace load
      else if (dict_component.type == 'load') {
        var explore = dict_component['plot_data'];
        explore.forEach(function(d){
            data_container.append("line")
              .attr('id', 'line_load')
              .attr("x1", 1000*d.data[0])
              .attr("y1", 1000*d.data[1])
              .attr("x2", 1000*d.data[2])
              .attr("y2", 1000*d.data[3])
              .attr("stroke-width", 1000*d.width)
              .attr("opacity", 0.7)
              .attr("stroke", 'blue')
              .attr("fill", "none")
              .attr("stroke-linecap", "round")
              .attr("marker-end", "url(#triangle_load)")
        })
      }

      // Trace contours
      else if (dict_component.type == 'contour') {
        var circles_data = [];
        var rectangles_data = [];
        var lines_data = [];
        var texts_data = [];
        var arcs_data = [];
        var wires_data = [];
        var area_data = [];

        var explore = dict_component['plot_data'];
        explore.forEach(function(d){
          if(d.type == "circle"){
            circles_data.push(d);
          }
          else if (d.type == "arc"){
            arcs_data.push(d);
            d.data.forEach(function(d){
              area_data.push([1000*d['x'],1000*d['y']])})
          }
          else if (d.type == "wire"){
            wires_data.push(d);
            d.data.forEach(function(d){
              area_data.push([1000*d['x'],1000*d['y']])})
          }
          else if (d.type == "line"){
            lines_data.push(d);
            area_data.push([1000*d.data[0],1000*d.data[1]]);
            area_data.push([1000*d.data[2],1000*d.data[3]]);
          }
          else if (d.type == "rect"){
            rectangles_data.push(d);
          }
        })
        arcs_data.forEach(function(d){
          var draw_data = []
          d.data.forEach(function(dd){
            draw_data.push([1000*dd['x'],1000*dd['y']])
          })
          var lineGenerator = d3.line();
          lineGenerator.curve(d3.curveCardinal.tension(0));
          var pathString = lineGenerator(draw_data);
          data_container.append("path")
              .attr('d', pathString)
              .attr('id', 'path_geom')
              .attr("stroke-width", d.size)
              .attr("stroke", d.color)
              .attr("fill", "none")
              .attr("stroke-width-init", d.size)
        })

        wires_data.forEach(function(d){
          var draw_data = []
          d.data.forEach(function(dd){
            draw_data.push([1000*dd['x'],1000*dd['y']])
          })
          var lineGenerator = d3.line();
          var pathString = lineGenerator(draw_data);
          data_container.append("path")
              .attr('d', pathString)
              .attr('id', 'path_geom')
              .attr("stroke-width", d.size)
              .attr("stroke", d.color)
              .attr("fill", "none")
              .attr("stroke-width-init", d.size)
        })

        circles_data.forEach(function(d){
          data_container.append("circle")
              .attr('r', 1000*d.r)
              .attr('cx', 1000*d.cx)
              .attr('cy', 1000*d.cy)
              .attr('id', 'circle_geom')
              .attr("stroke-width", d.size)
              .attr("stroke-width-init", d.size)
              .attr("stroke", 'black')
              .attr("fill", "none")
              .attr("pointer-events", 'all')
              //.on("mouseover", highlight)
              //.on("mouseout", deemphasize);
        })

        lines_data.forEach(function(d){
          var stroke_color = d3.rgb(d.color[0]*255, d.color[1]*255, d.color[2]*255);
          data_container.append("line")
            .attr('id', 'line_geom')
            .attr("x1", 1000*d.data[0])
            .attr("y1", 1000*d.data[1])
            .attr("x2", 1000*d.data[2])
            .attr("y2", 1000*d.data[3])
            .attr("stroke-width", d.size)
            .attr("stroke", d.color)
            // .attr("opacity", d.opacity)
            .attr("fill", "none")
            .attr("stroke-width-init", d.size)
            // .attr("stroke-linecap", "round");
        })

        rectangles_data.forEach(function(d){
          var stroke_color = d3.rgb(d.stroke[0]*255, d.stroke[1]*255, d.stroke[2]*255);
          data_container.append("rect")
          .attr('id', 'rect_geom')
          .attr("x", d.x*1000)
          .attr("y", d.y*1000)
          .attr("width", d.width*1000)
          .attr("height", d.height*1000)
          .attr("stroke", stroke_color.toString())
          .attr("stroke-width", d.size)
          .attr("stroke-width-init", d.size)
          .attr("fill", "none")
          .attr("pointer-events", "all")
          .attr("data", d);
          // .on("mouseover", highlight)
          // .on("mouseout", deemphasize);
        })

        // data_container.selectAll("rect")
        //     .data(rectangles_data)
        //     .enter()
        //     .append("rect")
        //     .attr('id', 'rect_geom')
        //     .attr("x", function(d){return d.x*1000;})
        //     .attr("y", function(d){return d.y*1000;})
        //     .attr("width", function(d){return d.width*1000;})
        //     .attr("height", function(d){return d.height*1000;})
        //     .attr("stroke", function(d){console.log(d);
        //       let color:any = d3.rgb(d.color[0]*256, d.color[1]*256, d.color[2]*256);
        //       return color;})
        //     .attr("fill", "none")
        //     .attr("pointer-events", "all")
        //     .on("mouseover", highlight)
        //     .on("mouseout", deemphasize);

        data_container.selectAll("text")
            .data(texts_data)
            .enter()
            .append("text")
            .attr("id", function(d){return d.text;})
            .attr("x", function(d){return d.x*1000;})
            .attr("y", function(d){return d.y*1000;})
            .attr("font-size", font_size)
            .attr('text-anchor', "middle")
            .text(function(d){return d.text;});

        if (dict_component['fill'] != 'none'){
          var areaGenerator = d3.area();
          var area = areaGenerator(area_data);
          data_container.append('path')
            .attr('id', 'essai')
            .attr('d', area)
            .attr("fill", "url(#diagonal-stripe-1)")
        }
      }
    })
      // if(d.type == "bspline"){ // Bspline type
      //   var lineGenerator = d3.line();
      //   lineGenerator.curve(d3.curveCardinal.tension(0));
      //   var pathString = lineGenerator(d.data);
      //
      //   if (d.group==1){
      //     d3.selectAll("#gear1").append('path')
      //       .attr('d', pathString)
      //       .attr("stroke-width", d.size)
      //       .attr("stroke", "blue")
      //       .attr("stroke-dasharray","none")
      //       .attr("class","gear")
      //       .attr("fill", "none");
      //   }
      //   else if (d.group==2){
      //     d3.selectAll("#gear2").append('path')
      //       .attr('d', pathString)
      //       .attr("stroke-width", d.size)
      //       .attr("stroke", "red")
      //       .attr("stroke-dasharray","none")
      //       .attr("class","gear")
      //       .attr("fill", "none");
      //   }
      // }
      // else if(d.type == "arrow"){ // Arrow type
        // var lineGenerator = d3.line();
        // var pathString = lineGenerator(d.data);
        // quoteComp+=1;
        // if (quoteComp % 2 ==0){
        //   var quoteOffset = "10%";
        // }
        // else {
        //   var quoteOffset = "70%";
        // }
        //
        // var quoteVal = Math.pow(Math.pow(d.data[0][0]-d.data[1][0],2)+Math.pow(d.data[0][1]-d.data[1][1],2),0.5);
        //
        // d3.selectAll("#Quote").append('path')
        //     .attr("d",pathString)
        //     .attr("stroke-width", d.size)
        //     .attr("stroke", "black")
        //     .attr("class","quote")
        //     .attr("fill", "none")
        //     .attr("marker-start", "url(#start)")
        //     .attr("marker-end", "url(#arrow)")
        //     .attr("id", "quote_"+quoteComp);
        //
        // d3.selectAll("#Quote")
        //     .append("text")
        //     .attr("class", "quoteText")
        //     .style('font-size', quoteFont)
        //     .attr("dy", -1)
        //     .append("textPath")
        //     .attr("xlink:href","#quote_"+quoteComp)
        //     .attr("startOffset",quoteOffset)
        //     .text("\u2300 "+precisionRound(quoteVal,2)+" mm");
      // }
      // else{
      //   console.log('Type '+d.type+' not implemented yet !');
      // }
    // });
    // }

    // Highlight function
    function highlight(d){
      d3.select(this)
        .attr("fill", "#e0e0e0")
        .attr("fill-opacity", 0.5)
        .attr("stroke", "#263238")
        .attr("stroke-width", this['stroke-width']*2/scale_factor); // Size à changer en stroke_width
    }

    // De-emphasize function
    function deemphasize(d){
      var stroke_color = d3.rgb(this.stroke[0]*255, this.stroke[1]*255, this.stroke[2]*255);
      d3.select(this)
        .attr("fill", "none")
        .attr("stroke", stroke_color.toString())
        .attr("stroke-width", this['stroke-width']/scale_factor); // Size à changer en stroke_width
    }
  }
}
