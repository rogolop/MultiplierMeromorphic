
AttachSpec("SingularitiesDim2/IntegralClosureDim2.spec");
Attach("MultiplierMeromorphic.m");

P<x, y> := LocalPolynomialRing(RationalField(), 2, "lglex");

// Bernstein-Sato...
// I := [P| y^3+x^5, x];
// I := [P| y^3+x^5, y];
// I := [P| y^3+x^5, 1];

// I := [P| x^5*y - y^3, x^8 + 2*x^5*y];

// I := [P| x^2-y^3, y^2-x^3];
// I := [P| x^3, y^2];
// I := [P| y^2, x^3];
// I := [P| x^3+y^2, y^2];
// I := [P| x^3+x*y^2, y^2];
// I := [P| (x^2-y^3)*(y^2-x^3), 1];
// I := [P| (x^2-y^3)*(y^2-x^3), x^2*y^2];

// I := [P| (x^2-y^3)*(x^2+y^3), y^2-x^3];
// I := [P| (x^2-y^3)^2, y^2-x^3];
// I := [P| (x^2-y^3)*x, 1];
// I := [P| (x^2-y^3)^2, 1];
// I := [P| (x^2-y^3)*y, 1];
// I := [P| y^2-x^3, (x^2-y^3)*(x^2+y^3)];
// I := [P| y*(x^2-y^3)*(x^2+y^3), y^2-x^3];
// I := [P| (x^9-y^4-x^3*y^3)^2 - x^16*y, y^2-x^3];

I := [P| (x^2-y^3)*(x^2+y^3), 1];
// I := [P| (x^2-y^3)^2, 1];
// I := [P| (x^2-y^3), 1];

// I := [P| (y^2-x^3)^5 + x^18, 1];
// I := [P| (y^2-x^3)^5 + x^18, y^2-x^3];
// I := [P| (y^2-x^3)^5 + x^18, (y^2-x^3)^2];
// I := [P| (y^2-x^3)^5 + x^18, (y^2-x^3)^3];
// I := [P| (y^2-x^3)^5 + x^18, (y^2-x^3)^4];
// I := [P| (y^2-x^3)^5 + x^18, (y^2-x^3)^5];

// I := [P| (y^2-x^3)^5 + x^18, y^2+x^3];

// I := [P| (y^2-x^3)^4 + x^8*y^5, 1];
// I := [P| (y^2-x^3)^4 + x^8*y^5, y^2-x^3];
// I := [P| (y^2-x^3)^4 + x^8*y^5, y^2+x^3];





f, g := Explode(I);
printf "f = %o\n", f;
printf "g = %o\n", g;
// printf "\nIdeal I=<%o, %o>\n", Basis(I)[1], Basis(I)[2];

// printf "\n### LogResolution\n";
// LogResolution(I: Coefficients:=true);

// printf "\n### LogResolutionMeromorphic\n";
// LogResolutionMeromorphic(I: Coefficients:=true);

// printf "\n### ProximityMatrix\n";
// ProximityMatrix(f: Coefficients := true, ExtraPoint := true);

// printf "\n### MultiplierIdeals I\n";
// S := MultiplierIdealsMeromorphic_I(I : MaxJN:=1);
// printf "%o\n", [t[2] : t in S];

// printf "\n### MultiplierIdeals f\n";
S := MultiplierIdealsMeromorphic(f, g : MinJN:=0, MaxJN:=1, ComputeIdeals:=false); /////////////////////
// printf "%o\n", [t[2] : t in S];
// S;


// f := (y^2-x^3)*x;
// f := (x^2-y^3)*(y^2-x^3);
// I := ideal<P | x^2*y^2, x^5, y^5, x*y^4, x^4*y>;

// printf "\n### MultiplierIdeals\n";
// printf "%o\n", [t[2] : t in S];
// S;
// S := MultiplierIdeals(I : MaxJN:=2);

// S := MultiplierIdeals(f : MaxJN:=3);
// S;
printf "\n%o\n", [t[2] : t in S]; ////////////////////////
den := LCM([Denominator(t[2]) : t in S]);
printf "\nNumerators:\n%o\n", [t[2]*den : t in S];

