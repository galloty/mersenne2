\\ Search for a primitive root whose order is 9*2^30
\\ such that primroot^{order/8) = 2^{(q - 1)/2} * (1, 1)

main()=
{
	q = 31; c = 9/2;
	\\q = 61; c = 9/2^31;
	p = 2^q - 1; 
	g = ffgen(Mod(1, p)*i^2 + Mod(1, p), 'i);
	r_8 = (1 + g) * 2^((q - 1)/2);
	r0 = ffprimroot(g)^((p - 1)/c);
	o = fforder(r0);

	j0 = 0; forstep (j = 1, 8, 2, if ((r0^j)^(o/8) == r_8, j0 = j;););

	d_min = p^2;
	r = r0^j0; r8 = r0^8;
	forstep (j = j0, o, 8,
		v = component(r, 1); d = v[2]^2 + v[3]^2;
		if (d < d_min,
			if (fforder(r) == o,
				d_min = d;
				print("primroot = ", r, ":");
				print("  j = ", j , " / ", o);
				print("  order = ", factor(o));
				print("  primroot^{order/2) = ", r^(o/2));
				print("  primroot^{order/4) = ", r^(o/4));
				print("  primroot^{order/8) = ", r^(o/8));
			);
		);

		r *= r8;
	);
}

main()
