#include <pari/pari.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <time.h>
#include "ServerBank.h"

using namespace std;

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
    pk1 = a;
    struct privatekey *b = new privatekey;
    b->p = p;
    b->q = q;
    b->vp = vp;
    b->vq = vq;
    sk1 = b;
}

GEN encrypt(GEN m, GEN n, GEN g, GEN h, GEN u, GEN l, GEN t){
    
    publickey *pk = new publickey;
    pk->n = n;
    pk->g = g;
    pk->h = h;
    pk->u = u;
    pk->l = l;
    pk->t = t;
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

GEN serverA(GEN x, GEN a, GEN n, GEN g, GEN h, GEN u, GEN l, GEN t, GEN p, GEN q, GEN vp, GEN vq, int ln){
    
    publickey *pk = new publickey;
    pk->n = n;
    pk->g = g;
    pk->h = h;
    pk->u = u;
    pk->l = l;
    pk->t = t;
    
    privatekey *sk = new privatekey;
    sk->p = p;
    sk->q = q;
    sk->vp = vp;
    sk->vq = vq;
    
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
		tmp = encrypt(stoi(c[i]), n, g, h, u, l, t);
		alpha.push_back(tmp);
	}
	gammap = serverB2Control(alpha, n, g, h, u, l, t, ln);
	for(int i=0; i<lengthm; i++){
		tmp = mulii(sk->vp, sk->vq);
  		tmp2 = Fp_pow(gammap[i], tmp, pk->n);
  		pari_printf("Decrypted text (should be 1 if plaintext = 0) = %Ps\n", tmp2);
	}

	return c;
}

