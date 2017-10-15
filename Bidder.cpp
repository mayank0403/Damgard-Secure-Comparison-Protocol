#include <pari/pari.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <time.h>
#include "ServerAuctionHouse.h"

using namespace std;

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
  int ln = lengthm;
  a = const_vec(lengthm, m);
  b = const_vec(lengthm, m);
  // Here (in the main function) we will consider the part of client.
  // tmp has the binary representation of the message.
  for(int i=1; i<=lengthm; i++){
  	tmp2 = randomi(pk1->u);
  	a[i] = itos(tmp2);
  	if(tmp[i]==0){
  		b[i] = itos(gsub(pk1->u, tmp2));
  	}
  	else{
  		tmp4 = gadd(pk1->u, tmp3);
  		b[i] = itos(gsub(tmp4, tmp2));
  	}
  }
  bool flg = false;
  for(int i=1; i<=lengthm; i++){
  	if((a[i]+b[i])%itos(pk1->u) != tmp[i]){
  		flg = true;
  		break;
  	}
  }
  
  //cout<<"\n--------------------------------------Transmitting to server----------------------------------------\n\n";

  tmp5 = serverB(x, b, pk1->n, pk1->g, pk1->h, pk1->u, pk1->l, pk1->t, ln);
  tmp4 = serverA(x, a, pk1->n, pk1->g, pk1->h, pk1->u, pk1->l, pk1->t, sk1->p, sk1->q, sk1->vp, sk1->vq, ln);
  
  pari_close();
  return 0;
}
