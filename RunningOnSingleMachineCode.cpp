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
} *pk;

struct privatekey{
	GEN p, q, vp, vq;
} *sk;

struct ll{
	GEN val;
	struct ll* next;
};

void keygen(){
  GEN t, l, u, vp, vq, tmp, tmp2, tbit, p, q, tmp3, n, h, g, tmp4;
  
  gp_read_stream(stdin);
  printf("l = "); l = gp_read_stream(stdin);
  printf("t(taken in decimal digits) = "); t = gp_read_stream(stdin);
  initial = clock();
  tmp2 = stoi(2);
  tmp = gadd(l, tmp2);
  u = gnextprime(tmp);
  tmp2 = stoi(10);
  // constructing 2 t-bit primes vp and vq such that vp|p-1 and vq|q-1
  tbit = powii(tmp2, t);
  tmp2 = stoi(2);
  vp = gnextprime(tbit);
  bool flag = false;
  tmp3 = stoi(1);
  for(int i=0; i<1000000; i++){
  	p = mulii(u, vp);
  	p = gadd(p, tmp3);
  	for(int j=0; j<=4; j++){
  		if(isprime(p)){
  			flag = true;
  			break;
  		}
  		p = gsub(p, tmp3);
  		p = mulii(p, tmp2);
  		p = gadd(p, tmp3);
  	}
  	if(flag){
  		break;
  	}
  	vp = gadd(vp, tmp3);
  	vp = gnextprime(vp);
  }
  tmp = gadd(vp, tmp3);
  vq = gnextprime(tmp);
  flag = false;
  for(int i=0; i<1000000; i++){
  	q = mulii(u, vq);
  	q = gadd(q, tmp3);
  	for(int j=0; j<=4; j++){
  		if(isprime(q)){
  			flag = true;
  			break;
  		}
  		q = gsub(q, tmp3);
  		q = mulii(q, tmp2);
  		q = gadd(q, tmp3);
  	}
  	if(flag){
  		break;
  	}
  	vq = gadd(vq, tmp3);
  	vq = gnextprime(vq);
  }
  n = mulii(p, q);
  cout<<"Number of bits in RSA modulus : "<<bit_prec(n)<<endl;
  tmp = gsub(p, tmp3);
  tmp = gsub(q, tmp3);
  // Selecting g and h
  // g and h should have value less than n and gcd with n = 1

  // Finding h first
  // order of h is vp.vq
  tmp2 = mulii(vp, vq);
  for (int i=20000000; i<50000000; i++){
  	tmp = stoi(i);
  	tmp4 = Fp_pow(tmp, tmp2, n);
  	if(equalii(tmp4, tmp3)){
  		break;
  	}
  }
  h = tmp;

  // Finding g secondly
  // order of g is u.vp.vq
  tmp2 = mulii(vp, vq);
  tmp2 = mulii(tmp2, u);
  for (int i=20000000; i<50000000; i++){
  	tmp = stoi(i);
  	tmp4 = Fp_pow(tmp, tmp2, n);
  	if(equalii(tmp4, tmp3)){
  		break;
  	}
  }
  g = tmp;
  struct publickey *a = new publickey;
  a->n = n;
  a->g = g;
  a->h = h;
  a->u = u;
  a->l = l;
  a->t = t;
  pk = a;
  struct privatekey *b = new privatekey;
  b->p = p;
  b->q = q;
  b->vp = vp;
  b->vq = vq;
  sk = b;
}

GEN encrypt(GEN m){
  GEN r, c, tmp, tmp2, tmp3, tmp4, tbit;
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
  r = tmp2;
  tmp = Fp_pow(pk->g, m, pk->n);
  tmp2 = Fp_pow(pk->h, r, pk->n);
  c = Fp_mul(tmp, tmp2, pk->n);
  return c;

}

GEN beta;
int cnt=0, ln;

GEN serverB2(GEN alpha){
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

GEN serverB(GEN x, GEN b){
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
	}
	beta = c;
	return c;
}


vector<GEN> serverB2Control(vector<GEN> alpha){
	vector<GEN> gammap;
	for(int i=0; i<ln; i++){
		cnt++;
		gammap.push_back(serverB2(alpha[i]));
	}
	random_shuffle ( gammap.begin(), gammap.end());
	return gammap;
}

GEN serverA(GEN x, GEN a){
	GEN w, c, tmp, tmp2, tmp3, tmp4, tmp5;
	tmp3 = stoi(1);
	tmp = stoi(0);
	int lengthm = lg(a)-1;
	w = const_vec(lengthm, stoi(1));
	c = const_vec(lengthm, stoi(1));
	for(int i=1; i<=lengthm; i++){
		tmp2 = gadd(stoi(x[i]), stoi(a[i]));
		tmp4 = mulii(stoi(x[i]), stoi(a[i]));
		tmp5 = stoi(2);
		tmp = mulii(tmp4, tmp5);
		w[i] = itos(gsub(tmp2, tmp));
	}
	for(int i=1; i<=lengthm; i++){
		tmp = stoi(0);
		for(int j=i+1; j<=lengthm; j++){
			tmp = gadd(tmp, stoi(w[j]));
		}
		tmp2 = gadd(tmp, stoi(x[i]));
		tmp2 = gadd(tmp2, tmp3);
		c[i] = itos(gsub(tmp2, stoi(a[i])));
		
	}
	vector<GEN> alpha, gammap;
	for(int i=1; i<=lengthm; i++){
		tmp = encrypt(stoi(c[i]));
		alpha.push_back(tmp);
	}
	gammap = serverB2Control(alpha);
	for(int i=0; i<lengthm; i++){
		tmp = mulii(sk->vp, sk->vq);
  		tmp2 = Fp_pow(gammap[i], tmp, pk->n);
  		pari_printf("Decrypted text (should be 1 if plaintext = 0) = %Ps\n", tmp2);
	}

	return c;
}

GEN dectobin(GEN a, int b){
	GEN ret = binary_zv(a);
	int len = lg(ret)-1;
	GEN ans, tmp;
	ans = const_vec(b, a);
	for(int i=b-len+1; i<=b; i++){
		tmp = stoi(ret[i-(b-len)]);
		ans[i] = itos(tmp);
	}
	for(int i=1; i<=b-len; i++){
		ans[i] = 0;
	}
	return ans;
}

int main()
{
  cout<<"NOTE : If the protocol gives incorrect result then increasing the size of l will make it work correctly\n\n";
  pari_init(500000000,2);
  cout<<"For 32 bit input use : type 1 and for 16 bit input type 2\n"<<endl;
  int in;
  cin>>in;
  keygen();
  afterkeygen = clock();
  GEN m, m1, x, x1, tmp, tmp2, tmp3, tmp4, tmp5, a, b;
  tmp3 = stoi(1);
  m = stoi(1);
  if(in == 1){
  	m1 = stoi(230000301);
  	x1 = stoi(230000000);
  	tmp = dectobin(m1, 32);
  	x = dectobin(x1, 32);	
  }
  else{
  	m1 = stoi(1200);
  	x1 = stoi(500);
  	tmp = dectobin(m1, 16);
  	x = dectobin(x1, 16);
  }
  tmp = vecsmall_reverse(tmp);
  x = vecsmall_reverse(x);
  int lengthm = lg(tmp)-1;
  ln = lengthm;
  a = const_vec(lengthm, m);
  b = const_vec(lengthm, m);
  // Here (in the main function) we will consider the part of client.
  // tmp has the binary representation of the message.
  for(int i=1; i<=lengthm; i++){
  	tmp2 = randomi(pk->u);
  	a[i] = itos(tmp2);
  	if(tmp[i]==0){
  		b[i] = itos(gsub(pk->u, tmp2));
  	}
  	else{
  		tmp4 = gadd(pk->u, tmp3);
  		b[i] = itos(gsub(tmp4, tmp2));
  	}
  }
  bool flg = false;
  for(int i=1; i<=lengthm; i++){
  	if((a[i]+b[i])%itos(pk->u) != tmp[i]){
  		flg = true;
  		break;
  	}
  }
  
  //cout<<"\n--------------------------------------Transmitting to server----------------------------------------\n\n";
  
  tmp5 = serverB(x, b);
  tmp4 = serverA(x, a);
  
  pari_close();
  return 0;
}
