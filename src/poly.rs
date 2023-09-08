use crate::{params::*, utils::*, ntt::*};

/* polynomial ring (make a struct?) arithm op modulo q, ntt domain */

// todo: i think only square matrices are needed here

#[derive(Copy, Clone)]
pub struct Poly
{
    pub coeff: [i16; N],
    pub ntt: bool,
}

impl Default for Poly
{
  fn default() -> Self
  {
    Poly { coeff: [0i16; N] , ntt: false}
  }
}

pub fn neg(a: &Poly) -> Poly
{
    let mut res = Poly::default();
    for i in 0..N 
    {
        res.coeff[i] = cmod(-a.coeff[i], Q as i16);
    }
    res.ntt = a.ntt;
    res    
}

pub fn add(a: &mut Poly, b: &Poly)
{
    assert_eq!(a.ntt, b.ntt);
    if a.ntt == b.ntt
    {
        for i in 0..N 
        {
            a.coeff[i] = cmod(a.coeff[i] + b.coeff[i], Q as i16);
        }
    }
}


pub fn ntt(a: &Poly) -> Poly
{
    let mut res = Poly::default();
    res.coeff = a.coeff;
    if ! a.ntt
    {
        ntt_(&mut res.coeff);
        res.ntt = true;
    }
    res
}

pub fn reduce(a: &mut Poly)// -> i32
{
    for i in 0..N 
    {
        a.coeff[i] = cmod(montgomery(a.coeff[i] as i32), Q as i16);
    }
}

pub fn inv_ntt(a: &Poly) -> Poly
{
    let mut res = Poly::default();
    res.coeff = a.coeff;
    if a.ntt
    {
        inv_ntt_(&mut res.coeff);
        reduce(&mut res);
        res.ntt = false;
    }
    res
}

/* scalar multiplication */
pub fn smult(s: i16, a: &Poly) -> Poly
{
    let mut res = Poly::default();
    for i in 0..N 
    {
        res.coeff[i] = s*a.coeff[i];
    }
    res.ntt = a.ntt;
    res
}

pub fn mult(a: &Poly, b: &Poly) -> Poly
{
    assert_eq!(a.ntt, b.ntt);
    assert_eq!(a.ntt, true);
    let mut res = Poly::default();
    for i in 0..N 
    {
        res.coeff[i] = montgomery((a.coeff[i] as i32) * b.coeff[i] as i32);
        let t = (((Q as i32 + res.coeff[i] as i32)<<16 ) % (Q as i32)) as i16; 
        res.coeff[i] = crate::utils::cmod(t, Q as i16); 
    }
    res.ntt = true;
    res
}

pub fn slow_mult(a: &Poly, b: &Poly) -> Poly
{
    assert_eq!(a.ntt, b.ntt);
    assert_eq!(a.ntt, false);
    let mut res = Poly::default();
    for i in 0..N 
    {
        for j in 0..N-i
        {
            res.coeff[i+j] += a.coeff[i] * b.coeff[j];           
        }
    } 
    for j in 1..N 
    {
        for i in N-j .. N 
        {
            res.coeff[i+j-N] -= a.coeff[i] * b.coeff[j];
        }
    }
    res.ntt = false;
    res
}


pub fn p_infnorm(p: &Poly) -> i16
{
    let mut norm = 0i16;
    for i in 0..N 
    {
        let t = p.coeff[i] - ((p.coeff[i] >> 15) & 2 * p.coeff[i]);
        if t >= norm 
        {
            norm = t;
        }
    }
    norm
}

/* vector of polynomials */
#[derive(Copy, Clone)]
pub struct VecPoly<const l: usize>
{
    pub poly: [Poly; l],
}


impl<const l:usize> Default for VecPoly<{l}>
{
  fn default() -> Self
  {
    VecPoly { poly: [Poly::default(); l] }
  }
}


pub fn inf_norm<const k: usize>(p: &VecPoly<{k}>) -> i16
{
    let mut norm = 0i16;
    for i in 0..k
    {
        if p_infnorm(&p.poly[i]) >= norm 
        {
            norm = p_infnorm(&p.poly[i]);
        }
    }
    norm
}

/* matrix of polynomials */
#[derive(Copy, Clone)]
pub struct Mat<const k: usize, const l:usize>
{
    pub vec: [VecPoly<l>; k],
}


impl<const k: usize, const l:usize> Default for Mat<{k}, {l}>
{
  fn default() -> Self
  {
    Mat { vec: [VecPoly::<l>::default(); k] }
  }
}

pub fn transp<const k: usize, const l: usize>(A: &Mat<{k}, {l}>) -> Mat<{l}, {k}>
{
    let mut res = Mat::<l, k>::default();
    for i in 0..l
    {
        for j in 0..k
        {
            res.vec[i].poly[j] = A.vec[j].poly[i];
        }
    }
    res
}

pub fn b_reduce(a: &mut Poly)
{
    for i in 0..N 
    {
        a.coeff[i] = barrett(a.coeff[i]);
    }
}

pub fn vb_reduce<const k: usize>(a: &mut VecPoly<{k}>) 
{
    for i in 0..k 
    {
        b_reduce(&mut a.poly[i]);
    }
}


pub fn Ntt<const k: usize>(a: &VecPoly<{k}>) -> VecPoly<{k}>
{
    let mut res = *a;//.copy();
    for i in 0..res.poly.len()
    {
        res.poly[i] = ntt(&a.poly[i]);
        res.poly[i].ntt = true;
    }
    res
}


pub fn inv_Ntt<const k: usize>(a: &VecPoly<{k}>) -> VecPoly<{k}>
{
    let mut res = VecPoly::<{k}>::default();
    for i in 0..res.poly.len()
    {
        res.poly[i] = inv_ntt(&a.poly[i]);
    }
    res
}


pub fn Neg<const k: usize>(a: &VecPoly<{k}>) -> VecPoly<{k}>
{
    let mut res = VecPoly::<{k}>::default();
    for i in 0..k 
    {
        res.poly[i] = neg(&a.poly[i]);
    }
    res    
}

pub fn sMult<const k: usize>(s: i16, a: &VecPoly<{k}>) -> VecPoly<{k}>
{
    let mut res = VecPoly::<{k}>::default();
    for i in 0..k 
    {
        res.poly[i] = smult(s, &a.poly[i]);
    }
    res  
}


pub fn v_add<const k: usize>(a: &VecPoly<k>, b: &VecPoly<k>) -> VecPoly<k>
{
    let mut res = VecPoly::<k>::default();
    for i in 0..k 
    {
        res.poly[i] = a.poly[i];
        add(&mut res.poly[i], &b.poly[i]);
    }
    res
}


pub fn v_sub<const k: usize>(a: &VecPoly<k>, b: &VecPoly<k>) -> VecPoly<k>
{
    let mut res = VecPoly::<k>::default();
    for i in 0..k 
    {
        res.poly[i] = a.poly[i];
        add(&mut res.poly[i], &smult(-1, &b.poly[i]));
    }
    res
}


pub fn p_mult_v<const k: usize>(p: &Poly, v: &VecPoly<{k}>) -> VecPoly<{k}>
{
    let mut res = VecPoly::<k>::default();
    for i in 0..k 
    {
        res.poly[i] = mult(&p, &v.poly[i]); 
    }
    Caddq(&res);

    res
}


/* scalar vector multiplication */
pub fn v_mult_v<const k:usize>(v: &VecPoly<{k}>, u: &VecPoly<{k}>) -> Poly
{
    let mut res = Poly::default();
    res.ntt = true;
    for i in 0..k 
    {
        add(&mut res, &mult(&v.poly[i], &u.poly[i])); 
    }
    res
}



/* matrix and vector multiplication */
pub fn m_mult_v<const k:usize, const l:usize>(A: &Mat<{k}, {l}>, s: &VecPoly<{l}>) -> VecPoly<{k}>
{
    let mut res = VecPoly::<{k}>::default();
    for i in 0..k 
    {
        res.poly[i] = v_mult_v(&A.vec[i], &s);
    }
    res
}

/* vector and matrix multiplication */
pub fn v_mult_m<const k:usize, const l:usize>(s: &VecPoly<{k}>, A: &Mat<{k}, {l}>) -> VecPoly<{l}>
{
    let mut res = VecPoly::<{l}>::default();
    for i in 0..l 
    {
        res.poly[i] = v_mult_v(&s, &transp(A).vec[i]);
    }
    res
}



/* utilities */
pub fn c_addq(a: &Poly) -> Poly
{
    let mut res = Poly::default();
    for i in 0..N 
    {
        res.coeff[i] = cmod(a.coeff[i], Q as i16);
    }
    res.ntt = a.ntt;
    res
}

pub fn _modq(a: &Poly) -> Poly
{
    let mut res = Poly::default();
    for i in 0..N 
    {
        res.coeff[i] = mod_a(a.coeff[i], Q as i16);
        
    }
    res.ntt = a.ntt;
    res
}


pub fn Caddq<const k: usize>(a: &VecPoly<k>) -> VecPoly<k>
{
    let mut res = VecPoly::<k>::default();
    for i in 0..k 
    {
        res.poly[i] = c_addq(&a.poly[i]);
    }
    res
}

pub fn modq<const k: usize>(a: &VecPoly<k>) -> VecPoly<k>
{
    let mut res = VecPoly::<k>::default();
    for i in 0..k 
    {
        res.poly[i] = _modq(&a.poly[i]);
    }
    res
}


pub fn compress(r: &Poly, d: i16) -> Poly
{
    let mut t = Poly::default();
    for i in 0..r.coeff.len()
    {
        t.coeff[i] = _compress(r.coeff[i], d);
    }
    t.ntt = r.ntt;
    t
}


pub fn decompress(r: &Poly, d: i16) -> Poly
{
    let mut t = Poly::default();
    for i in 0..r.coeff.len()
    {
        t.coeff[i] = _decompress(r.coeff[i], d);
    }
    t.ntt = r.ntt;
    t
}

pub fn Compress(r: &VecPoly<K>, d: i16) -> VecPoly<K>
{
    let mut t = VecPoly::<K>::default();
    for i in 0..r.poly.len()
    {
        t.poly[i] = compress(&r.poly[i], d);
    }
    t
}


pub fn Decompress(r: &VecPoly<K>, d: i16) -> VecPoly<K>
{
    let mut t = VecPoly::<K>::default();
    for i in 0..r.poly.len()
    {
        t.poly[i] = decompress(&r.poly[i], d);
    }
    t
}