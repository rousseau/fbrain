struct RBF_fn {
	virtual Doub rbf(Doub r) = 0;
};

struct RBF2_fn {
  virtual Doub rbf(const Doub *p1, const Doub *p2) = 0;
};

struct RBF_interp {
	Int dim, n;
	const MatDoub &pts;
	const VecDoub &vals;
	VecDoub w;
	RBF2_fn &fn;
	Bool norm;

//	RBF_interp(MatDoub_I &ptss, VecDoub_I &valss, RBF_fn &func, Bool nrbf=false)
//	: dim(ptss.ncols()), n(ptss.nrows()) , pts(ptss), vals(valss),
//	w(n), fn(func), norm(nrbf) {
//		Int i,j;
//		Doub sum;
//		MatDoub rbf(n,n);
//		VecDoub rhs(n);
//		for (i=0;i<n;i++) {
//			sum = 0.;
//			for (j=0;j<n;j++) {
//				sum += (rbf[i][j] = fn.rbf(rad(&pts[i][0],&pts[j][0])));
//			}
//			if (norm) rhs[i] = sum*vals[i];
//			else rhs[i] = vals[i];
//		}
//		LUdcmp lu(rbf);
//		lu.solve(rhs,w);
//	}

  RBF_interp(MatDoub_I &ptss, VecDoub_I &valss, RBF2_fn &func, Bool nrbf=false)
  : dim(ptss.ncols()), n(ptss.nrows()) , pts(ptss), vals(valss),
  w(n), fn(func), norm(nrbf) {
    Int i,j;
    Doub sum;
    MatDoub rbf(n,n);
    VecDoub rhs(n);

//    std::cout << " Puntos = " << std::endl;
//    for (i=0;i<n;i++) {
//      for (j=0;j<dim;j++) {
//        std::cout << pts[i][j] << " " ;
//      }
//      std::cout << std::endl;
//    }
//
//    std::cout << " Valores = " << std::endl;
//    for (i=0;i<n;i++) {
//      std::cout << vals[i] << " ";
//    }
//    std::cout << std::endl;

//    std::cout << " Matriz phi = " << std::endl;
    for (i=0;i<n;i++) {
      sum = 0.;
      for (j=0;j<n;j++) {
        sum += (rbf[i][j] = fn.rbf(&pts[i][0],&pts[j][0]));
//        std::cout << rbf[i][j] << " " ;
      }
//      std::cout << std::endl;
      if (norm) rhs[i] = sum*vals[i];
      else rhs[i] = vals[i];
    }
    LUdcmp lu(rbf);
    lu.solve(rhs,w);
  }

//	Doub interp(VecDoub_I &pt) {
//		Doub fval, sum=0., sumw=0.;
//		if (pt.size() != dim) throw("RBF_interp bad pt size");
//		for (Int i=0;i<n;i++) {
//			fval = fn.rbf(rad(&pt[0],&pts[i][0]));
//			sumw += w[i]*fval;
//			sum += fval;
//		}
//		return norm ? sumw/sum : sumw;
//	}

  Doub interp(VecDoub_I &pt) {
    Doub fval, sum=0., sumw=0.;
    if (pt.size() != dim) throw("RBF_interp bad pt size");
    for (Int i=0;i<n;i++) {
      fval = fn.rbf(&pt[0],&pts[i][0]);
      sumw += w[i]*fval;
      sum += fval;
//      std::cout <<  w[i] << " ";
    }
 //   std::cout <<  std::endl;
    return norm ? sumw/sum : sumw;
  }

	Doub rad(const Doub *p1, const Doub *p2) {
		Doub sum = 0.;
		for (Int i=0;i<dim;i++) sum += SQR(p1[i]-p2[i]);
		return sqrt(sum);
	}
};
struct RBF_multiquadric : RBF_fn {
	Doub r02;
	RBF_multiquadric(Doub scale=1.) : r02(SQR(scale)) {}
	Doub rbf(Doub r) { return sqrt(SQR(r)+r02); }
};

struct RBF_thinplate : RBF_fn {
	Doub r0;
	RBF_thinplate(Doub scale=1.) : r0(scale) {}
	Doub rbf(Doub r) { return r <= 0. ? 0. : SQR(r)*log(r/r0); }
};

struct RBF_gauss : RBF_fn {
	Doub r0;
	RBF_gauss(Doub scale=1.) : r0(scale) {}
	Doub rbf(Doub r) { return exp(-0.5*SQR(r/r0)); }
};

struct RBF2_gauss : RBF2_fn {
  Doub r0_spa;
  Doub r0_ang;
  RBF2_gauss(Doub scale1=1, Doub scale2=1) : r0_spa(scale1),r0_ang(scale2) {}
  Doub rbf(const Doub *p1, const Doub *p2)
  {
    Doub r_spa = 0;
    Doub r_ang = 0;

    double q1[3];
    double q2[3];

    for (Int i=0;i<3;i++) r_spa += SQR(p1[i]-p2[i]);
    r_spa = sqrt( r_spa );

    q1[0] = std::sin(p1[3]) * std::cos(p1[4]);//m_GradientTableCartesian[1][0];
    q1[1] = std::sin(p1[3]) * std::sin(p1[4]);
    q1[2] = std::cos(p1[3]);

    q2[0] = std::sin(p2[3]) * std::cos(p2[4]);//m_GradientTableCartesian[1][0];
    q2[1] = std::sin(p2[3]) * std::sin(p2[4]);
    q2[2] = std::cos(p2[3]);

    for (Int i=0;i<3;i++) r_ang += q1[i]*q2[i];
    // value can be higher than 1 at machine precision
    if ( abs(r_ang) > 1  )
    {
//      std::cout << " entro !!! : " << r_ang << " " << abs(r_ang) << " " << (abs(r_ang) - 1) << std::endl;
      r_ang = 1;
    }
    r_ang = std::acos(abs( r_ang ));

    return exp(-0.5*SQR(r_spa/r0_spa))*exp(-0.5*SQR(r_ang/r0_ang));
  }
};

struct RBF_inversemultiquadric : RBF_fn {
	Doub r02;
	RBF_inversemultiquadric(Doub scale=1.) : r02(SQR(scale)) {}
	Doub rbf(Doub r) { return 1./sqrt(SQR(r)+r02); }
};
struct Shep_interp {
	Int dim, n;
	const MatDoub &pts;
	const VecDoub &vals;
	Doub pneg;

	Shep_interp(MatDoub_I &ptss, VecDoub_I &valss, Doub p=2.)
	: dim(ptss.ncols()), n(ptss.nrows()) , pts(ptss),
	vals(valss), pneg(-p) {}

	Doub interp(VecDoub_I &pt) {
		Doub r, w, sum=0., sumw=0.;
		if (pt.size() != dim) throw("RBF_interp bad pt size");
		for (Int i=0;i<n;i++) {
			if ((r=rad(&pt[0],&pts[i][0])) == 0.) return vals[i];
			sum += (w = pow(r,pneg));
			sumw += w*vals[i];
		}
		return sumw/sum;
	}

	Doub rad(const Doub *p1, const Doub *p2) {
		Doub sum = 0.;
		for (Int i=0;i<dim;i++) sum += SQR(p1[i]-p2[i]);
		return sqrt(sum);
	}
};
