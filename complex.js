/* complex.js
 * A super-simple complex number type and a few related functions.
 */

var eps = function()   // Machine epsilon utility function
{
  var x = 1;
  var y = x;
  while(1 + x > 1)
  {
    y = x;
    x = x/2;
  }
  return y;
};

/* A complex type.
 * x: Either a complex value or a real value.
 * y: Either a real value representing the imaginary part, or undefined.
 *
 * If x is a complex then the new complex value is a copy of x and
 * the argument y is ignored.
 * If x is a number then the it's set to the real part of the new value.
 * If y is undefined and x is real, then the new imaginary part will be
 * set to zero.
 */
Complex = function(x, y)
{
  if(x.re == undefined)
  {
    this.re = x;
    if(y == undefined)
    {
      this.im = 0;
    } else
    {
      this.im = y;
    }
  } else
  {
    this.re = x.re;
    this.im = x.im;
  }
};

Complex.prototype.add = function(z)  // Sum
{
  if(z.re == undefined) z = new Complex(z,0);
  return new Complex(this.re + z.re, this.im + z.im);
};

Complex.prototype.sub = function(z)  // Difference
{
  if(z.re == undefined) z = new Complex(z,0);
  return new Complex(this.re - z.re, this.im - z.im);
};

Complex.prototype.mul = function(z)  // Product
{
  if(z.re == undefined) z = new Complex(z,0);
  return new Complex(this.re * z.re - this.im * z.im,
                     this.re * z.im + this.im * z.re);
};

Complex.prototype.div = function(z)  // Quotient
{
  if(z.re == undefined) z = new Complex(z,0);
  return this.mul(z.conj()).mul(1/Math.pow(z.mod(),2));
};

Complex.prototype.conj = function()  // Conjugate
{
  return new Complex(this.re, -this.im);
};

Complex.prototype.mod = function()  // Modulus
{
  return Math.sqrt(this.re*this.re + this.im*this.im);
};

Complex.prototype.arg = function()  // Argument
{
  return Math.atan2(this.im, this.re);
};

Complex.prototype.pow = function(q)  // Exponentiation
{
  var r = Math.pow(this.mod(),q);
  var t = this.arg();
  return new Complex(r*Math.cos(t*q), r*Math.sin(t*q));
};

Complex.prototype.sign = function()
{
  if(this.mod() < 2*eps())
  {
    return 1;
  }
  return this.div(this.mod());
};


/* Miscellaneous complex linear algebra routines...
 *
 * The following basic routines assume that matrices are stored as flat arrays
 * of column-major data.  Use these carefully, the eigenvalue routines in
 * particular are not well-optimized. See the in-line notes.
 */

/* convenience function for vector inner product y'x */
var ip = function(x,y)
{
  var n = x.length;
  if(y.length !=n) throw(new Error("the vectors are not the same length"));
  return mm(t(y,n,1),x,1,n,1)[0];
};

/* vector axpy */
var axpy = function(a, x, y)
{
  if(y == undefined)
  {
    var y = zeros(x.length);
  }
  for(var j=0;j<x.length;++j)
  {
    y[j] = y[j].add(x[j].mul(a));
  }
  return y;
};

/* (Conjugate) transpose of an m x n matrix or vector A.
 * Set conjugate=false to just return the transpose, otherwise returns
 * the conjugate transpose.
 */
t = function(A, m, n, conjugate)
{
  if(conjugate == undefined)
  {
    conjugate = true;
  }
  var B = new Array(m*n);
  for(var i=0;i<m;++i)
  {
    for(var j=0;j<n;++j)
    {
      if(conjugate)
      {
        B[i*n + j] = A[i + j*m].conj();
      } else {
        B[i*n + j] = A[i + j*m];
      }
    }
  }
  return B;
};

/* Matrix multiply.
 * Let A and B be vectors of complex values. Interpret A as an m*k matrix
 * whose values are stored in column-major order, and similarly interpret
 * B as a k*n matrix. Form the m*n complex-valued matrix product AB.
 */
mm = function(A, B, m, k, n)
{
  if((A.length != m*k) || (B.length != k*n))
  {
    throw(new Error("Non-conforming lengths"));
  }
  var h,i,j;
  var C = zeros(m*n);
  for(h=0;h<n;++h)
  {
    for(i=0;i<m;++i)
    {
      for(j=0;j<k;++j)
      {
        C[h*m + i] = C[h*m+i].add(A[j*m+i].mul(B[h*k+j]));
      }
    }
  }
  return C;
};

/* Assume that A is a square matrix presented in column-major order.
 * Compute its trace.
 */
trace = function(A)
{
  var n = Math.sqrt(A.length);
  var z = new Complex(0,0);
  for(var j=0;j<n;++j)
  {
    z = z.add(A[j + j*n]);
  }
  return z;
};

/* Compute the eigenvalues of a 2x2 matrix directly via the quadratic
 * formula applied to its characteristic polynomial.
 */
var eigs2 = function(A)
{
  var B = A[0].add(A[3]);
  var C = A[0].mul(A[3]).sub(A[2].mul(A[1]));
  var d = (B.mul(B).sub(C.mul(4))).pow(0.5);
  var l = [];
  B = B.mul(-1);
  l.push(B.add(d).div(2));
  l.push(B.sub(d).div(2));
  return l;
};

/* Compute the eigenvalues of a 3x3 matrix directly by solving for the
 * roots of its characteristic polynomial.
 */
var eigs3 = function(A)
{
  if(Math.sqrt(A.length) !=3) throw(new Error("eigs3 requires a 3x3 matrix"));
  var d = A[6].mul(A[4]).mul(A[2]).
          add(A[3].mul(A[1]).mul(A[8])).
          add(A[0].mul(A[7]).mul(A[5])).
          sub(A[0].mul(A[4]).mul(A[8])).
          sub(A[3].mul(A[7]).mul(A[2])).
          sub(A[6].mul(A[1]).mul(A[5]));
  var b = trace(A).mul(-1);
  var c = trace(mm(A,A,3,3,3)).sub(b.mul(b)).mul(-0.5);
  var a = new Complex(1, 0);
  var d0 = b.mul(b).sub(a.mul(c).mul(3));
  var d1 = b.mul(b).mul(b).mul(2).
           sub(a.mul(b).mul(c).mul(9)).
           add(a.mul(a).mul(d).mul(27));
  var C = d1.mul(d1).sub(d0.mul(d0).mul(d0).mul(4)).pow(1/2).add(d1).mul(1/2).pow(1/3);
  var z = new Complex(-1/3,0);
  var u1 = new Complex(1,0);
  var u2 = new Complex(-1/2, Math.sqrt(3)/2);
  var u3 = new Complex(-1/2, -Math.sqrt(3)/2);
  u1 = u1.mul(C);
  u1 = d0.div(u1).add(u1).add(b).mul(z);
  u2 = u2.mul(C);
  u2 = d0.div(u2).add(u2).add(b).mul(z);
  u3 = u3.mul(C);
  u3 = d0.div(u3).add(u3).add(b).mul(z);
  return [u1,u2,u3];
};

/* Return a zero array */
var zeros = function(n)
{
  return Array.prototype.map.call(Array.apply(null,new Array(n)),
           function(){return new Complex(0,0);});
};

/* Return an nxn identity matrix */
var identity = function(n)
{
  var I = new Array(n*n);
  for(var k=0;k<n*n;++k)
  {
    I[k] = new Complex(0,0);
  }
  for(var k=0;k<n;++k)
  {
    I[k*n+k] = new Complex(1,0);
  }
  return I;
};

var show = function(A,m,n)
{
  for(var i=0;i<m;i++)
  {
    var line = "";
    for(var j=0;j<n;j++)
    {
      var k = j*m + i;
      line = line + A[k].re.toFixed(2) + " + " + A[k].im.toFixed(2) + "i\t";
    }
    console.log(line);
  }
};

/* Givens rotation for zero-indexed element in row i > column j of square
 * matrix A.  A must be column-major ordered. See lawn 148. Return a rotation
 * matrix of the same size as A.
 */
var givens = function(A,i,j)
{
  var n = Math.sqrt(A.length);
  var a = A[j*n + j];
  var b = A[j*n + i]; // The element to be rotated to zero
// TODO Handle exceptional cases.
  var r = Math.sqrt(a.mod()*a.mod() + b.mod()*b.mod());
  var c = new Complex(a.mod()/r, 0);
  var s = b.conj().div(r).mul(a.sign());
  var G = identity(n);
  G[j*n + j] = c;
  G[j*n + i] = s.conj().mul(-1);
  G[i*n + j] = s;
  G[i*n + i] = c;
  return G;
};

/* QR decomposition of a general matrix A. Returns two matrices Q and R
 * so that A = QR.
 */
var qr = function(A)
{
  var n = Math.sqrt(A.length);
  var epsilon = 2*eps();
  var Q = identity(n);
  var R = A;
  for(var j=0;j<n-1;++j)
  {
    for(var i=j+1;i<n;++i)
    {
      if(R[j*n+i].mod()>epsilon) // The entry to rotate to zero
      {
        G = givens(R,i,j);
        R = mm(G,R,n,n,n);
        Q = mm(Q,t(G,n,n),n,n,n);
      }
    }
  }
  var ans = new Object();
  ans["Q"] = Q;
  ans["R"] = R;
  return ans;
};

/* QR decomposition of an upper Hessenberg matrix A. Returns two matrices Q and
 * R so that A = QR and Q is unitary and R is upper triangular.
 */
var qrhess = function(A)
{
  var n = Math.sqrt(A.length);
  var epsilon = 2*eps();
  var Q = identity(n);
  var R = A;
  var i;
  for(var j=0;j<n-1;++j)
  {
    i = j+1;
    if(R[j*n+i].mod()>epsilon) // The entry to rotate to zero
    {
      G = givens(R,i,j);
      R = mm(G,R,n,n,n);
      Q = mm(Q,t(G,n,n),n,n,n);
    }
  }
  var ans = new Object();
  ans["Q"] = Q;
  ans["R"] = R;
  return ans;
};

var norm2 = function(v)
{
  return Math.sqrt(ip(v,v).re);
};

/* Arnoldi Hessenberg reduction of a square matrix A.
 * Returns matrices V, H where:
 * V is unitary
 * H is upper Hessenberg
 * AV = VH.
 */
var arnoldihess = function(A)
{
  var n = Math.sqrt(A.length);
  var V = zeros(n*(n+1));
  var H = zeros(n*n);
  var v = zeros(n);
  var w, wn;
  for(var k=0;k<n;++k)
  {
    v[k] = new Complex(Math.random()-0.5,0);
  }
  v = axpy(1/norm2(v), v);
  for(var k=0;k<n;++k)
  {
    V[k] = v[k];
  }
  var k;
  var j = 0;
  while(j<n)
  {
    w = mm(A, V.slice(j*n,j*n+n),n,n,1);
    for(k=0;k<=j;++k)
    {
      H[j*n+k] = ip(w,V.slice(k*n,k*n+n));
      w = axpy(H[j*n+k].mul(-1), V.slice(k*n,k*n+n), w);
    }
    wn = norm2(w);
    H[j*n + j + 1] = new Complex(wn,0);
    w = axpy(1/wn,w);
    j = j + 1;
    for(k=0;k<n;++k)
    {
      V[j*n + k] = w[k];
    }
  }
  ans = new Object();
  ans["V"] = V.slice(0,n*n);
  ans["H"] = H.slice(0,n*n);
  return ans;
};

// Roots of a3*z^3 + a2*z^2 + a1*z1 + a0
var cubic = function(a3,a2,a1,a0)
{
  var A3 = new Complex(a3);
  var A2 = new Complex(a2);
  var A1 = new Complex(a1);
  var A0 = new Complex(a0);
  var d0 = A2.pow(2).sub(A3.mul(A1).mul(3));
  var d1 = A2.pow(3).mul(2)
             .sub(A1.mul(A2).mul(A3).mul(9))
             .add(A3.pow(2).mul(A0).mul(27));
  var C  = (d1.add((d1.pow(2).sub(d0.pow(3).mul(4))).pow(1/2)).mul(1/2)).pow(1/3);
  var z = A3.pow(-1).mul(-1/3);
  var u1 = new Complex(1,0);
  var u2 = new Complex(-1/2, Math.sqrt(3)/2);
  var u3 = new Complex(-1/2, -Math.sqrt(3)/2);
  u1 = u1.mul(C);
  u1 = d0.div(u1).add(u1).add(A2).mul(z);
  u2 = u2.mul(C);
  u2 = d0.div(u2).add(u2).add(A2).mul(z);
  u3 = u3.mul(C);
  u3 = d0.div(u3).add(u3).add(A2).mul(z);
  return [u1,u2,u3];
};

// cubic resolvent roots, assuming real-valued coefficients
var resolvent = function(a3, a2, a1, a0)
{
  var A3 = 1;
  var A2 = -a2;
  var A1 = a1*a3 - 4*a0;
  var A0 = 4*a2*a0 - a1*a1 -a3*a3*a0;
  var ans = cubic(A3,A2,A1,A0);
  return ans;
};

/* Roots of the monic fourth order polynomial
 * z^4 + a3*z^3 + a2*z^2 + a1*z1 + a0
 * with real-valued coefficients.
 */
var quartic = function(a3,a2,a1,a0)
{
  var tol = 0.000001;
  var tol2 = 0.00001;
  var rz = resolvent(a3,a2,a1,a0);
  var y1 = rz.filter(function(z) {return(Math.abs(z.im)<tol);})[0];
  var D, E;
  var R = y1.add(0.25*a3*a3).sub(a2).pow(0.5);
  if(R.mod() < tol2)  // Tolerance here is a problem XXX
  {
    D = y1.pow(2).sub(4*a0).pow(0.5).mul(2).sub(2*a2).add(0.75*a3*a3).pow(0.5);
    E = y1.pow(2).sub(4*a0).pow(0.5).mul(-2).sub(2*a2).add(0.75*a3*a3).pow(0.5);
  } else
  {
    D = R.pow(2).mul(-1).add(0.75*a3*a3 - 2*a2).
                 add(R.pow(-1).mul(0.25*(4*a3*a2 - 8*a1 - a3*a3*a3))).pow(0.5);
    E = R.pow(2).mul(-1).add(0.75*a3*a3 - 2*a2).
                 sub(R.pow(-1).mul(0.25*(4*a3*a2 - 8*a1 - a3*a3*a3))).pow(0.5);
  }
  var z1 = R.mul(0.5).sub(a3/4).add(D.mul(0.5));
  var z2 = R.mul(0.5).sub(a3/4).sub(D.mul(0.5));
  var z3 = R.mul(-0.5).sub(a3/4).add(E.mul(0.5));
  var z4 = R.mul(-0.5).sub(a3/4).sub(E.mul(0.5));
  return [z1, z2, z3, z4];
};
