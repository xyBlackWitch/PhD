
void test_poca_neutral_neutral(){

  // This is the exact(!) skewed line POCA
   //
   TVector3 vertex;
   TVector3 av;
   TVector3 ap;
   TVector3 bv;
   TVector3 bp;

   vertex.SetXYZ(0.,0.,0.);
   bv.SetXYZ(-7.,2.,-3.);
   bp.SetXYZ(0.,1.,2.);
   av.SetXYZ(-3.,-3.,3.);
   ap.SetXYZ(1.,2.,1.);


   TVector3   u = ap.Unit();
   TVector3   v = bp.Unit();
   TVector3   w = av - bv;
   double     a = u.Mag(); //dot(u,u);  == 1       // always >= 0
   double     b = u.Dot(v); //dot(u,v); Skalarprodukt
   double     c = v.Mag(); //dot(v,v);  == 1       // always >= 0
   double     d = u.Dot(w); //dot(u,w); Skalarprodukt
   double     e = v.Dot(w); //dot(v,w);
   double     D = a*c - b*b;        // always >= 0
   double     sc, tc;

   cout << "D: " << D << " a: " << a << " b: " << b <<" c: " << c << " d: " << d << " e: " << e << endl;


  // compute the line parameters of the two closest points
   if (D < 1e-9) {          // the lines are almost parallel
	 sc = 0.0;
	 tc = e;    // use the largest denominator
   }
   else {
       sc = (b*e - c*d) / D; //determination of s in g: x1 = av + s*ap
       tc = (a*e - b*d) / D; //determination of t in h: x2 = bv + t*bp
       	  	  	  	  	  	//using (x1-x2)*u = 0 and (x1-x2)*v = 0
   }
   cout << " sc: " << sc << " tc: " << tc << endl;

   // get the difference of the two closest points
   //TVector3   dP = w + (sc * u) - (tc * v);  // =  L1(sc) - L2(tc)
   TVector3   pocaa = av + sc*u;
   TVector3   pocab = bv + tc*v;
   vertex = 0.5*(pocaa+pocab);
   TVector3   diff = pocaa-pocab;

   cout << "Diff " << diff.Mag() << endl;   // return the closest distance
   cout << "Vertex Position : (" << vertex.X() << ", " << vertex.Y() << ", " << vertex.Z() << ") " << endl;
   cout << "Pocaa Position : (" << pocaa.X() << ", " << pocaa.Y() << ", " << pocaa.Z() << ") " << endl;
   cout << "Pocab Position : (" << pocab.X() << ", " << pocab.Y() << ", " << pocab.Z() << ") " << endl;

}

