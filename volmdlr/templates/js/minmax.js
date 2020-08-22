var evaluate_minmax = function(data, coeff_pixel) {
    this.evaluate = function(){
      var minX = 0;
      var maxX = 0;
      var minY = 0;
      var maxY = 0;
      data.forEach(function(dict_component){
        var texts_data = [];
        if (dict_component['type'] == 'line'){
          if(dict_component['data'][0]*coeff_pixel < minX){
            minX = dict_component['data'][0]*coeff_pixel;
          }
          if(dict_component['data'][1]*coeff_pixel < minY){
            minY = dict_component['data'][1]*coeff_pixel;
          }
          if(dict_component['data'][2]*coeff_pixel > maxX){
            maxX = dict_component['data'][2]*coeff_pixel;
          }
          if(dict_component['data'][3]*coeff_pixel > maxY){
            maxY = dict_component['data'][3]*coeff_pixel;
          }
        }
        if (dict_component['type'] == 'contour'){
        var explore = dict_component['plot_data_primitives'];
        explore.forEach(function(d){
          if (d['type'] == "circle" || d['type'] == "arc"){
            if((d['cx'] - d['r'])*coeff_pixel < minX){
              minX = (d['cx'] - d['r'])*coeff_pixel;
            }
            if((d['cx'] + d['r'])*coeff_pixel > maxX){
              maxX = (d['cx'] + d['r'])*coeff_pixel;
            }
            if((d['cy'] - d['r'])*coeff_pixel < minY){
              minY = (d['cy'] - d['r'])*coeff_pixel;
            }
            if((d['cy'] + d['r'])*coeff_pixel > maxY){
              maxY = (d['cy'] + d['r'])*coeff_pixel;
            }
          }
          else if(d['type'] == "rect"){
            if(d['x']*coeff_pixel < minX){
              minX = d['x']*coeff_pixel;
            }
            if((d['x'] + d['width'])*coeff_pixel > maxX){
              maxX = (d['x'] + d['width'])*coeff_pixel;
            }
            if(d['y']*coeff_pixel < minY){
              minY = d['y']*coeff_pixel;
            }
            if((d['y'] + d['height'])*coeff_pixel > maxY){
              maxY = (d['y'] + d['height'])*coeff_pixel;
            }
          }
        })}
      })
      return [minX, maxX, minY, maxY]
    }

  var sol = this.evaluate()
  this.minX = sol[0]
  this.maxX = sol[1]
  this.minY = sol[2]
  this.maxY = sol[3]

}
