#include <pari/pari.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <time.h>

using namespace std;

clock_t initial, final, afterkeygen;

GEN
extgcd(GEN A, GEN B, GEN *U, GEN *V)
{
  pari_sp av = avma;
  GEN ux = gen_1, vx = gen_0, a = A, b = B;

  if (typ(a) != t_INT) pari_err_TYPE("extgcd",a);
  if (typ(b) != t_INT) pari_err_TYPE("extgcd",b);
  if (signe(a) < 0) { a = negi(a); ux = negi(ux); }
  while (!gequal0(b))
  {
    GEN r, q = dvmdii(a, b, &r), v = vx;

    vx = subii(ux, mulii(q, vx));
    ux = v; a = b; b = r;
  }
  *U = ux;
  *V = diviiexact( subii(a, mulii(A,ux)), B );
  gerepileall(av, 3, &a, U, V); return a;
}

// create a function to generate keys.

struct publickey{
	GEN n, g, h, u, l, t;
} *pk1;

struct privatekey{
	GEN p, q, vp, vq;
} *sk1;

struct ll{
	GEN val;
	struct ll* next;
};

GEN beta;
int cnt=0;

GEN serverB2(GEN alpha, GEN n, GEN g, GEN h, GEN u, GEN l, GEN t, int ln){
    
    publickey *pk = new publickey;
    pk->n = n;
    pk->g = g;
    pk->h = h;
    pk->u = u;
    pk->l = l;
    pk->t = t;
	GEN tmp, tmp2, tmp3, tmp4, s, sdash, tbit, e;
	tmp2 = stoi(10);
  	// constructing 2 t-bit primes vp and vq such that vp|p-1 and vq|q-1
  	tbit = powii(tmp2, pk->t);
  	tmp3 = stoi(1);
  	tmp2 = stoi(2);
  	tmp4 = mulii(pk->t, tmp2);
  	tmp4 = gadd(tmp4, tmp3);
  	tmp2 = stoi(10);
  	tbit = powii(tmp2, tmp4);
  	tmp2 = randomi(tbit);
  	tmp2 = gadd(tmp2, tbit);
  	sdash = tmp2;
  	s = randomi(pk->u);
  	while(itos(s)==0){
  		s = randomi(pk->u);
  	}
  	tmp2 = Fp_pow(pk->g, stoi(beta[cnt]), pk->n);
  	tmp = Fp_mul(alpha, tmp2, pk->n);
  	tmp2 = Fp_pow(tmp, s, pk->n);
  	tmp4 = Fp_pow(pk->h, sdash, pk->n);
  	e = Fp_mul(tmp2, tmp4, pk->n);
  	return e;
}

GEN serverB(GEN x, GEN b, GEN n, GEN g, GEN h, GEN u, GEN l, GEN t, int ln){
    publickey *pk = new publickey;
    pk->n = n;
    pk->g = g;
    pk->h = h;
    pk->u = u;
    pk->l = l;
    pk->t = t;
	GEN w, c, tmp, tmp2, tmp4, tmp5;
	tmp = stoi(0);
	int lengthm = lg(b)-1;
	w = const_vec(lengthm, stoi(1));
	c = const_vec(lengthm, stoi(1));
	for(int i=1; i<=lengthm; i++){
		tmp2 = gadd(stoi(0), stoi(b[i]));
		tmp4 = mulii(stoi(x[i]), stoi(b[i]));
		tmp5 = stoi(2);
		tmp = mulii(tmp4, tmp5);
		w[i] = itos(gsub(tmp2, tmp));
	}
	for(int i=1; i<=lengthm; i++){
		tmp = stoi(0);
		for(int j=i+1; j<=lengthm; j++){
			tmp = gadd(tmp, stoi(w[j]));
		}
		tmp2 = gadd(tmp, stoi(0));
		tmp2 = gadd(tmp2, stoi(0));
		c[i] = itos(gsub(tmp2, stoi(b[i])));
        cout<<c[i]<<" ";
	}
    cout<<endl;
	beta = c;
    
	return c;
}


vector<GEN> serverB2Control(vector<GEN> alpha, GEN n, GEN g, GEN h, GEN u, GEN l, GEN t, int ln){
	vector<GEN> gammap;
	for(int i=0; i<ln; i++){
		cnt++;
		gammap.push_back(serverB2(alpha[i], n, g, h, u, l, t, ln));
	}
	random_shuffle ( gammap.begin(), gammap.end());
	return gammap;
}

