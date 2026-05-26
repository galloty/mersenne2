\\ Search for a primitive root whose order is 9*2^30
\\ such that primroot^{order/8) = 2^{(q - 1)/2} * (1, 1)

main()=
{
	q = 61; c = 9/2^31;
	\\q = 31; c = 9/2;
	p = 2^q - 1; 
	g = ffgen(Mod(1, p)*i^2 + Mod(1, p), 'i);
	d_min = p^2;
	for (k = 1, oo,
		pr = ffprimroot(g);
		if (pr^((p^2 - 1)/8) == (1 + g) * 2^((q - 1)/2),
			r = pr^((p - 1)/c);
			v = component(r, 1); d = v[2]^2 + v[3]^2;
			if (d < d_min,
				d_min = d;
				o = fforder(r);
				print("primroot = ", r, ":");
				print("  order = ", factor(o));
				print("  primroot^{order/4) = ", r^(o/4));
				print("  primroot^{order/8) = ", r^(o/8));
			);
		);
	);
}

main()
