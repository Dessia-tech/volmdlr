
var hatching = function(max_size, nb_hatch, lineWidth){
  this.max_size = max_size
  this.nb_hatch = nb_hatch
  this.lineWidth = lineWidth

  this.generate = function(){
    var p_hatch = document.createElement("canvas");
    p_hatch.width = this.max_size;
    p_hatch.height = this.max_size;
    var pctx = p_hatch.getContext("2d");
    pctx.lineCap = 'square';
    pctx.beginPath();
    var pos_x = - Math.pow(Math.pow(this.max_size,2)/2, 0.5)
    var pos_y = Math.pow(Math.pow(this.max_size,2)/2, 0.5)
    for (var i = 0; i <= 2*this.nb_hatch; i++) {
      pos_x = pos_x + this.max_size/this.nb_hatch
      pos_y = pos_y - this.max_size/this.nb_hatch
      pctx.moveTo(pos_x, pos_y);
      pctx.lineTo(pos_x + this.max_size, pos_y + this.max_size);
    }
    pctx.stroke();
    return p_hatch
  }
}
