/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: Elisabeth
 *
 * Created on 4 avril 2017, 03:19
 */

#include <cstdlib>
#include <iostream>

using namespace std;

void aff(unsigned n)
{
	unsigned bit = 0 ;
	unsigned mask = 1 ;
	for (int i = 0 ; i < 32 ; ++i)
	{
		bit = (n & mask) >> i ;
		printf("%d", bit) ;
		mask <<= 1 ;
	}
}

/* Un polynôme s'écrira sous la forme : a0 + a1X + a2 X² + etc ???*/
/* Polynome[0]=a0, Polynome[1]=a1, ... */

/*uint8_t* Polynome; 
int degre;*/

// rev de p
/*void rev_table(uint8_t &p, int d){
    int i;
    uint8_t k;
    
    for(i=0; i<d+1; i++){
        k=p[i];
        p[i]=p[d-i];
        p[d-i]=k;
    }
}*/

// d= deg
unsigned int rev(unsigned int p, unsigned int d){
   int r=0;
   int i, k;
   
   for(i=0; i<=d; i++){
       k=(p&1);
       k=k<<(d-i);
       r=r|k;
       p>>=1;
       
       /*cout << "k=" << k << " r=" << r << endl;
       cout << "k="; aff(k); cout << endl;
       cout << " r="; aff(r); cout << endl;
       cout << "p="; aff(p); cout << endl;*/
   }
   
   return r;
}

unsigned int BezoutBinaire(unsigned int *pu, unsigned int *pv, unsigned int a, unsigned int b){
    int u, v, r, s, p, q;
    int d=0;
    
    while(((a&1)==0)&&((b&1)==0)){
        a>>=1; 
        b>>=1;
        d++;
    }
    
    u=1; r=0;
    v=0; s=1;
    p=a;
    q=b;
    
    while(q!=0){
        while((p&1)==0){
            p>>=1;
            
            if(((u&1)==0) && ((v&1)==0)) {
                u>>=1;
                v>>=1;
            }
            
            else{
                u=(u-b)>>1; v=(v+a)>>1;
            }
            
            while ((q&1)==0) {
                q>>=1;
                if(((r&1)==0) && ((s&1)==0)) {
                    s>>=1;
                    r>>=1;
                }
                else {
                    r=(r-b)>>1;
                    s=(s+a)>>1;
                }
            }
            
            if(p>q) {
                p-=q;
                v-=s;
                u-=r;
            }
            else {
                q-=p;
                s-=v;
                r-=u;
            }
        }
        
        *pu=u;
        *pv=v;
        
        return p<<d;
    }
}

unsigned int Bezout(unsigned int *pu, unsigned int *pv, unsigned int a, unsigned int b){
    int d=0, u=1, v=0, r=0, s=1, p=a, q=b;
    
    while((a&1==0) &&  ((b&1)==0)){
        a>>=1;
        b>>=1;
        d++;
    }
    
    while(q!=0){
        while((p&1)==0){
            p>>=1;
            if(((u&1)==0) && ((v&1)==0)){
                u>>=1;
                v>>=1;
            }
            else{
                u=(u-b)>>1;
                v=(v+a)>>1;
            }
        }
        
        while((q&1)==0){
            q>>=1;
            if(((r&1)==0) && ((s&1)==0)) 
                {s>>=1; r>>=1;}
            else{
                r=(r-b)>>1;
                s=(s+a)>>1;
            }
        }
        
        if(p>q){
            p-=q;
            v-=s;
            u-=r;
        }
        else{
            q-=p;
            s-=v;
            r-=u;
        }
    }
    
    // affectation des résultats et retour de la fonction
    *pu=u;
    *pv=v;
    
    return p<<d;
}

// algo 1 : calcul du polynôme minimal d'une suite
unsigned int pol_min(unsigned int d, unsigned int u){
    unsigned int s, t;
    unsigned int x=1<<(2*d); // x^{2d}
    
    Bezout(&t, &s, u, x);
    
    return rev(t, d);
}




/*
 
 * voir algo trouvé
 
 */

// Affiche un unsigned int en binaire :
void print_b(unsigned int a){
    int l=sizeof(unsigned int)*8-1; // position du dernier bit d'un unsigned int
    
    while(l!=0 && ((a>>l)&1)==0) l--; // on n'affiche pas les bits à 0 au début
    
    // On affiche les bits 1 à 1;
    while(l!=0){
        cout << ((a>>l)&1);
        l--;
    }
    cout << (a&1); // on affiche le dernier bit;
    
    cout << endl;
}

int mulmod(int a, int n, int m){
    return (a*n)% m;
}

int chiffre(int a, int i){
    int k=0;
    while(k++<i) a=a/10;
    return a%10;
}

int mettre_chiffre(int &a, int b, int i){
    /*int k=0, debut=a, fin;
    while(k++ < i) debut = debut/10;
    fin = a - debut;
    
    debut = debut - debut/10 + b;
    
    k=0;
    while(k++ < i) debut = debut*10;
   
    debut=debut+fin;
    return debut;*/
    
    int k=0, m=1;
    
    while(k++ < i) m*=10;
    
    a=a-((a/m)%10)*m+b*m;
    
    return a;
}

int puissance_de_10(int a){
    int k=0, tmp;
    
    if(a==0) return -1;
    
    while((tmp=a/10)!=0) {a=tmp; k++;}
    
    return k;
}

int diff_polynomes(int a, int b, int p){
    int res=0;
    int i, k;

    //cout << "puiss :" << puissance_de_10(a) << endl;
    
    for(i=puissance_de_10(a); i>=0; i--){
        //cout << "a=" << chiffre(a, i) << endl;
        //cout << "b=" << chiffre(b, i) << endl;
        k=(chiffre(a, i)+p-chiffre(b, i))%p;
        //cout << "k=" << k << endl;
        mettre_chiffre(res, k, i);
    }
    
    return res;
}

int mult_polynomes(int a, int b, int p){
    int ka=puissance_de_10(a);
    int kb=puissance_de_10(b);
    int res=0;
    int i, j;
    
    for(i=ka; i>=0; i--){
        for(j=kb; j>=0; j--){
            //cout << chiffre(a,i) << " " << chiffre(b,j) << endl;
            mettre_chiffre(res, (chiffre(res, i+j)+mulmod(chiffre(a, i), chiffre(b, j), p))%p, i+j);
        }
    }
    
    return res;    
}

int division(int a, int b, int p/*, int *reste=NULL*/){
    int kb=puissance_de_10(b);
    int q=0, r=a, kr=puissance_de_10(a);
    int i;
    int tmp=0;
    int stop =0;
    
    
    while(kr>=kb && stop++<10){
        for(i=0; i<p; i++){
            if(mulmod(i, chiffre(b, kb), p) == chiffre(r, kr)){
                tmp=0;
                mettre_chiffre(tmp, i, kr-kb);
                mettre_chiffre(q, i, kr-kb);
                /*cout << "r=" << r << endl;
                cout << "tmp=" << tmp << endl;*/
                r=diff_polynomes(r, mult_polynomes(b, tmp, p), p);
                /*cout << "on retire : " << mult_polynomes(b, tmp, p);
                cout << "r=" << r << endl;*/
                kr=puissance_de_10(r);
                
                break;
            }
        }
    }    
    
    /**reste=r;*/
    
    return q;
}

// p=2, 3, 5 ou 7 seulement :
int bezout(int *u, int *v, int a, int b, int p){
    int u0=1, v0=0, u1=0, v1=1, ut, vt;
    int r0=a, r1=b, rt;
    int k;
    
    //cout << "1 - dans bezout :" << u0 << " " << v0 << endl;
    
    while(r1!=0){
        k=division(r0, r1, p);
        
        cout << "k: " << k << endl;
        
        // colonne des restes :
        rt=diff_polynomes(r0, mult_polynomes(k, r1, p), p);
        r0=r1; 
        r1=rt;

        // colonne des u : 
        ut=diff_polynomes(u0, mult_polynomes(k, u1, p), p);
        u0=u1;
        u1=ut;

        // colonne des v :
        vt=diff_polynomes(v0, mult_polynomes(k, v1, p), p);
        v0=v1;
        v1=vt;
        
        cout << "dans bezout 0 :" << u0 << " " << v0 << endl;
        cout << "r0=" << r0 << endl;
        cout << "dans bezout 1 :" << u1 << " " << v1 << endl;
        cout << "r1=" << r1 << endl;
    }
    
    *u=u0;
    cout << "dans bezout 0 :" << u0 << " " << v0 << endl;
    cout << "dans bezout 1 :" << u1 << " " << v1 << endl;  
    
    *v=v0;
    
    return r0;
}

/*int bezout(uint8_t **u, uint8_t **v, uint8_t *a, uint8_t *b){
    uint8_t *u0, *v0, *u1, *v1;
    
    
}*/

int main(int argc, char** argv) {
    //unsigned int polynome=100010; 
    //unsigned int p2=0b110110;
    //cout << p2 << endl; 
    
    //cout << rev(p2, 6) << endl;
    
    //aff(pol_min(3, p2));

    int u, v;
    
    cout << "bezout : " << bezout(&u, &v, 1000000, 32403, 5) << endl;
    cout << "u=" << u << " v=" << v << endl;
    
    //cout << "mult = " << mult_polynomes(32403, 221, 5) << endl<<endl;
    //cout << "div : " << division(1000000, 32403, 5) << endl;
    //cout << "diff = " << diff_polynomes(1000000, 1430100, 5) << endl;
    //cout << "div : " << division(42, 3, 5) << endl;
    
    return 0;
}

