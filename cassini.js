/*
 * Copyright (c) 2014, B. W. Lewis All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 * 
 * 3. Neither the name of the copyright holder nor the names of its
 * contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * 
 */

// XXX Consider replacing maps with for loops.
// Apparently, imperative-style programming in JS is better optimized than
// declarative. See http://jsperf.com/map-vs-native-for-loop/7

var eigenvalues_showing = false;
var gershgorin_showing = false;
var cassini_showing = false;

/* Complex number generator */
var gencomplex = function()
{
  var x = Math.floor(Math.random()*10 - 5);
  var y = Math.floor(Math.random()*10 - 5);
  var v = x + " ";
  if(y<0)
  {
    v = v + " - " + Math.abs(y) + "i";
  } else if(y>0)
  {
    v = v + " + " + y + "i";
  }
  return v;
};

/* Parse a string representation of a complex number, returning
 * a Complex typed value.
XXX Improve this to recognize single i value
 */
var parse_complex = function(s)
{
  var v = s.replace(/([0-9]+)-/,"$1 -")
           .replace(/([0-9]+)\+/,"$1 +")
           .replace(/- +([0-9.]+)/,"-$1")
           .replace(/\+ +([0-9.]+)/,"+$1")
           .replace("i","").replace("I","")
           .replace(/^ +/,"").split(/ +/);
  var x = parseFloat(v[0]);
  var y = 0;
  if(v.length > 1)
  {
    y = parseFloat(v[1]);
    if(isNaN(y)) y = 0;
  }
  return new Complex(x,y);
};

/* Return a vector of complex numbers from the UI matrix. Values are returned
 * in column-major order.
 */
var get_values = function ()
{
  return t(d3.selectAll(".matrix")[0].map(function(x){return parse_complex(x.value);}),3,3,false);
};

var off_diag_row_sums = function (A)
{
  var n = Math.sqrt(A.length);
  var ans = [];
  for(var j=0;j<n;++j)
  {
    var s = 0;
    for(var k=0;k<n;++k)
    {
      if(k!=j)
      {
        s = s + A[k*n + j].mod();
      }
    }
    ans.push(s);
  }
  return ans;
};

var cls = function()
{
// Remove existing plot elements
  d3.selectAll(".labels").remove();
  d3.selectAll(".gershgorin").remove();
  d3.selectAll(".eigenvalues").remove();
  d3.selectAll(".cassini").remove();
  eigenvalues_showing = false;
  gershgorin_showing = false;
  cassini_showing = false;
};

var gersh = function(main, xax, yax, data, DELAY)
{
// Draw the Gershgorin discs
  var ggn = main.append("svg:g")
                .selectAll("gershgorin")
                .data(data)
                .enter().append("svg:circle")
                .attr("cx", function (d) { return xax(d[0]); } )
                .attr("cy", function (d) { return yax(d[1]); } )
                .attr("r", function(d) {return xax(d[2]) - xax(0);})
                .attr("class","gershgorin")
                .style("opacity","0")
                .style("fill", function(d){return d[3];})
                .transition().style("opacity","1")
                .duration(200).delay(DELAY);
  var label2 = main.append("svg:g")
                   .append("svg:text")
                   .attr("x","275").attr("y","0")
                   .attr("text-anchor","end")
                   .attr("class","labels")
                   .attr("cursor","pointer")
                   .style("font-family","sans-serif")
                   .style("stroke","#55a")
                   .style("opacity","0")
                   .text("Gershgorin discs")
                   .on('click', function()
                     {
                       if(gershgorin_showing)
                       {
                         d3.selectAll(".gershgorin").remove();
                         gershgorin_showing = false;
                       } else {
                         draw(500,2,0);
                       }
                     })
                   .transition().style("opacity","1")
                   .duration(200).delay(DELAY);
  gershgorin_showing = true;
}

var draw_eigs = function(main, xax, yax, A)
{
// Draw the eigenvalues
  var lambda = eigs3(A);
  main.append("svg:g")
      .selectAll("eigenvalues")
      .data(lambda)
      .enter().append("svg:line")
      .attr("x1", function (d) { return xax(d.re) - 5; } )
      .attr("x2", function (d) { return xax(d.re) + 5; } )
      .attr("y1", function (d) { return yax(d.im); } )
      .attr("y2", function (d) { return yax(d.im); } )
      .attr("class","eigenvalues")
      .style("stroke-width", "1")
      .style("stroke", "#000");
  main.append("svg:g")
      .selectAll("eigenvalues")
      .data(lambda)
      .enter().append("svg:line")
      .attr("x1", function (d) { return xax(d.re); } )
      .attr("x2", function (d) { return xax(d.re); } )
      .attr("y1", function (d) { return yax(d.im) - 5; } )
      .attr("y2", function (d) { return yax(d.im) + 5; } )
      .attr("class","eigenvalues")
      .style("stroke-width", "1")
      .style("stroke", "#000");
  var label1 = main.append("svg:g")
                   .append("svg:text")
                   .attr("x","100").attr("y","0")
                   .attr("text-anchor","end")
                   .attr("class","labels")
                   .attr("cursor","pointer")
                   .style("stroke","#555")
                   .style("opacity","1")
                   .style("font-family","sans-serif")
                   .on('click', function()
                     {
                       if(eigenvalues_showing)
                       {
                         d3.selectAll(".eigenvalues").remove();
                         eigenvalues_showing = false;
                       } else {
                         draw(500,1,0);
                       }
                     })
                   .text("Eigenvalues");
  eigenvalues_showing = true;
}

// Brauer's ovals of Cassini...
var draw_cassini = function(main, xax, yax, A, N, DELAY)
{
  var line = d3.svg.line()
               .x(function(d) { return xax(d.re); })
               .y(function(d) { return yax(d.im); })
               .interpolate("linear");

  var interp = function(A)
    {
      var o;
      for(j=0;j<A.length;++j)
      {
        o = main.append("svg:g").append("svg:path")
                .attr("d", line(A[j])).attr("class","cassini")
                .style("stroke-width","1")
                .style("stroke","#a00")
                .style("fill","#a00")
                .style("opacity","0");
        o.transition().style("opacity","0.2").duration(200).delay(DELAY);
      }
    };
  interp(brauer(A,0,1,N));
  interp(brauer(A,0,2,N));
  interp(brauer(A,1,2,N));

  var label3 = main.append("svg:g")
                   .append("svg:text")
                   .attr("x","500").attr("y","0")
                   .attr("text-anchor","end")
                   .attr("cursor","pointer")
                   .style("stroke","#a55")
                   .attr("class","labels")
                   .style("font-family","sans-serif")
                   .style("opacity","0")
                   .text("Brauer's ovals of Cassini")
                   .on('click', function()
                     {
                       if(cassini_showing)
                       {
                         d3.selectAll(".cassini").remove();
                         cassini_showing = false;
                       } else {
                         draw(500,3,0);
                       }
                     })
                   .transition().style("opacity","1")
                   .duration(200).delay(DELAY);
  cassini_showing = true;
}

var draw = function(N, which, DELAY)
{
  var A = get_values();
  var pad = 1;
  var radii = off_diag_row_sums(A);
  var data = [[A[0].re, A[0].im, radii[0], "rgba(0,0,255,0.5)"],
              [A[4].re, A[4].im, radii[1], "rgba(50,255,0,0.5)"],
              [A[8].re, A[8].im, radii[2], "rgba(255,255,0,0.5)"]];
  var lim = data.map(function(v){return v[0] - v[2] - pad;}).
                 concat(data.map(function(v){return v[0] + v[2] + pad;})).
                 concat(data.map(function(v){return v[1] - v[2] - pad;})).
                 concat(data.map(function(v){return v[1] + v[2] + pad;}));
  var pmin = d3.min(lim);
  var pmax = d3.max(lim);
// Reset the axes (width, height are global variables)
  var xax = d3.scale.linear().domain([pmin,pmax]).range([0,width]);
  var yax = d3.scale.linear().domain([pmin,pmax]).range([height,0]);
  if(which == 0)
  {
    cls();
    d3.select('.yaxis').call(d3.svg.axis().scale(yax).orient("left"));
    d3.select('.xaxis').call(d3.svg.axis().scale(xax).orient("bottom"));
  }

  if(which % 2 == 0) gersh(main,xax,yax, data, DELAY);
  if(which % 3 == 0) draw_cassini(main,xax,yax, A, N, 2*DELAY);
  if(which % 1 == 0) draw_eigs(main,xax,yax,A);

};

var rot = function(theta)
{
  var R = new Array(9);
  for(var k=0;k<9;++k)
  {
    R[k] = new Complex(0,0);
  }
  for(var k=0;k<3;++k)
  {
    R[k*3+k] = new Complex(Math.cos(theta),Math.sin(theta));
  }
  return R;
};

var shift = function (X, z)
{
  var n = Math.sqrt(X.length);
  var Y = X;
  for(var k=0;k<n;++k)
  {
    Y[k*n + k] = X[k*n + k].sub(z);
  }
  return Y;
};

var brauer = function(A,i,j,N)
{
  var theta = Array.apply(null,new Array(N)).map(
                     function(x,i) {return 2*Math.PI*i/N;});
  var r = off_diag_row_sums(A);
  var tol = Math.sqrt(eps());

  var z = A[i*3 + i];  // Diagonal element A[i,i]
  var r1 = -Math.atan(z.im/z.re);
  if(isNaN(r1))
  {
    r1 = 0;
  }
  var W1 = rot(r1);
  var A1 = mm(W1, A, 3, 3, 3);  // Rotate A[i,i] to the real line
  var s1 = A1[i*3 + i];
  var A2 = shift(A1,s1); // Shift A[i,i] to zero

  var w = A2[j*3 + j];  // Diagonal element A2[j,j]
  var r2 = -Math.atan(w.im/w.re);
  if(isNaN(r2))
  {
    r2 = 0;
  }
  var W2 = rot(r2);
  var A3 = mm(W2, A2, 3, 3, 3); // Rotate A[j,j] to the real line
  var s2 = A3[j*3 + j].re/2; // Diagonal element A3[j,j] divided by 2
  var A4 = shift(A3,s2);
  var p = r[i]*r[j];

  var solve = function(t)
       { 
         var b = real_root(t,s2,p);
         return b.map(function(w)
               {
                 var z = new Complex(w*Math.cos(t), w*Math.sin(t));
                 var flag = z.re>0;
                 var theta = z.add(s2).arg();
                 var gamma = z.sub(s2).arg();
// undo the rotations and shifts...
                 z = W1[0].conj().mul((W2[0].conj().mul(z.add(s2))).add(s1));
                 z.flag = flag;
                 z.theta = theta;
                 z.gamma = gamma;
                 return z;
               });
       };

  var u = theta.map(solve);
// There are two possible cases:
// Two roots define a single oval.
// Four roots defines two ovals.
  var m = d3.max(u.map(function(x) {return x.length;}));
  if(m<=2)
  {
    var X0 = [].concat.apply([],u).filter(function(z) {return !isNaN(z.re) && !isNaN(z.im)});
    X = [X0.sort(function(a,b){return a.arg() - b.arg();})];
    return X;
  } else
  {
    var X = [
[].concat.apply([],u.filter(function(v){return v.length==4;})).filter(function(z) {return !isNaN(z.re) && !isNaN(z.im) && z.flag}).sort(function(a,b){return a.theta - b.theta;}),
[].concat.apply([],u.filter(function(v){return v.length==4;})).filter(function(z) {return !isNaN(z.re) && !isNaN(z.im) && !z.flag}).sort(function(a,b){return a.gamma - b.gamma;})
];
    return X;
  }
};

/* Find the real-valued root(s) of the monic polynomial
 * x^4 - 2*a^2*cos(2*t)*x^2 + a^4 - K^2,
 * where 0 <= t < 2pi, a,K are real-valued and K > 0.
 */
real_root = function(t,a,K)
{
  var tol = 0.0001;
  var a2 = -2*a*a*Math.cos(2*t);
  var a0 = a*a*a*a - K*K;
  var l = quartic(0.0,a2,0.0,a0);
  var ans = l.filter(function(z)
            {
              return (Math.abs(z.im)<tol); 
            }).map(function(w){return w.re;});
  return ans;
};
