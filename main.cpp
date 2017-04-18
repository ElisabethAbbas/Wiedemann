#include <cstdlib>
#include <iostream>

using namespace std;

// (a+b) mod p
int addmod(int a, int b, int p){
    return (a+b)%p;
}

// (a*b) mod p
int mulmod(int a, int b, int p){
    return (a*b)%p;
}

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

// récupère de le coefficient de la puissance 10^i-ème dans a
int get(int a, int i){
    int k=0;
    while(k++<i) a=a/10;
    return a%10;
}

// remplace le coefficient de la puissance 10^i-ème par b dans a
int set(int &a, int b, int i){
    int k=0, m=1;
    
    while(k++ < i) m*=10;
    
    a=a-((a/m)%10)*m+b*m;
    
    return a;
}

// renvoie le degré du polynôme de a, 
// pour le polynôme nul, on renverra -1
int degre(int a){
    int k=0, tmp;
    
    if(a==0) return -1;
    while((tmp=a/10)!=0) {a=tmp; k++;}
    
    return k;
}

// renvoie le polynôme "reversal" pour un corps Z/pZ avec p = 2, 3, 5 ou 7
int rev(int p, int d){
   int i, r=0;
   
   for(i=0; i<=d; i++)
        set(r, get(p, i), d-i);
   
   return r;
}

int add_p(int a, int b, int p){
    int res=0;
    int i, k, d_a, d_b;
   
    // on récupère le degré de a et de b :
    d_a=degre(a);
    d_b=degre(b);
    
    // on met dans i le max de {degré de a; degré de b}
    if (d_a>d_b) 
        i=d_a;
    else 
        i=d_b;
    
    // on soustrait les coefficients de a et de b et on les place dans res :
    for(; i>=0; i--){
        k=(get(a, i)+get(b, i))%p;
        set(res, k, i);
    }
    
    return res;
}

// renvoie la différence du polynôme a par le polynôme b, dans le corps Z/pZ, avec p= 2, 3, 5 ou 7
int diff_p(int a, int b, int p){
    int res=0;
    int i, k, d_a, d_b;
   
    // on récupère le degré de a et de b :
    d_a=degre(a);
    d_b=degre(b);
    
    // on met dans i le max de {degré de a; degré de b}
    if (d_a>d_b) 
        i=d_a;
    else 
        i=d_b;
    
    // on soustrait les coefficients de a et de b et on les place dans res :
    for(; i>=0; i--){
        k=(get(a, i)+p-get(b, i))%p;
        set(res, k, i);
    }
    
    return res;
}

// renvoie la multiplication du polynôme a par le polynôme b, dans le corps Z/pZ, avec p= 2, 3, 5 ou 7
int mult_p(int a, int b, int p){
    int ka=degre(a);
    int kb=degre(b);
    int res=0;
    int i, j;
    
    // on multiplie les coefficients et on les places dans res :
    for(i=ka; i>=0; i--)
        for(j=kb; j>=0; j--){
            set(res, (get(res, i+j)+mulmod(get(a, i), get(b, j), p))%p, i+j);
            //cout << "i+j=" << i+j << " get=" << get(res, i+j) << endl;
        }
    
    return res;    
}

// division du polynôme a par le polynôme b, dans le corps Z/pZ, avec p= 2, 3, 5 ou 7
int div_p(int a, int b, int p){
    int kb=degre(b);
    int q=0, r=a, kr=degre(a);
    int i;
    int tmp=0;    
    
    //cout << "a=" << a << " b=" << b << endl;
    
    // tant que le degré du dividende est plus grand ou égal à celui du diviseur
    while(kr>=kb){
        // on cherche le coefficient qui, multiplié au diviseur, donne celui du dividende :
        for(i=0; i<p; i++){
            if(mulmod(i, get(b, kb), p) == get(r, kr)){
                // on met à jour une variable tempioraire :
                tmp=0;
                set(tmp, i, kr-kb);
                // on met à jour le quotient : 
                set(q, i, kr-kb);
                
                // on soustrait (comme dans l'algorithme) pour obtenir le nouveau dividende :
                r=diff_p(r, mult_p(b, tmp, p), p);
                // on met à jour le degré du dividende :
                kr=degre(r);
                break;
            }
        }
    }    
    
    return q;
}

// renvoie le pgcd et met à jour les coefficients u et v, 
// sur les polynômes a et b et dans un corps Z/pZ, avec p=2, 3, 5 ou 7
int bezout(int *u, int *v, int a, int b, int p){
    int u0=1, v0=0, u1=0, v1=1, ut, vt;
    int r0=a, r1=b, rt;
    int q;
    
    while(r1!=0){
        q=div_p(r0, r1, p);
        
        //cout << "q = " << q << endl;
        
        // colonne des restes :
        rt=diff_p(r0, mult_p(q, r1, p), p);
        r0=r1; 
        r1=rt;
        
        // colonne des u : 
        //cout << "u0 = " << u0 << endl;
        //cout << "u1 = " << u1 << endl;
        ut=diff_p(u0, mult_p(q, u1, p), p);
        u0=u1;
        u1=ut;

        // colonne des v :
        //cout << "v0 = " << v0 << endl;
        //cout << "v1 = " << v1 << endl;
        vt=diff_p(v0, mult_p(q, v1, p), p);
        v0=v1;
        v1=vt;
    }
    
    // mise à jour de u et v
    *u=u0;    
    *v=v0;
    
    return r0;
}

// renvoie le pgcd et met à jour les coefficients u et v, 
// sur les polynômes a et b et dans un corps Z/pZ, avec p=2, 3, 5 ou 7
// avec les conditions de l'algo 1
int bezout_algo1(int *u, int *v, int a, int b, int p){
    int u0=1, v0=0, u1=0, v1=1, ut, vt;
    int r0=a, r1=b, rt;
    int q;
    int d_n = degre(b)/2;
    
    while(r1!=0){
        q=div_p(r0, r1, p);
        
        // colonne des restes :
        rt=diff_p(r0, mult_p(q, r1, p), p);
        r0=r1; 
        r1=rt;
        
        // colonne des u : 
        ut=diff_p(u0, mult_p(q, u1, p), p);
        u0=u1;
        u1=ut;

        // colonne des v :
        vt=diff_p(v0, mult_p(q, v1, p), p);
        v0=v1;
        v1=vt;
        
        // mise à jour de u et v
        if(degre(u1)<=d_n && degre(r1)<d_n){
            *u=u1;    
            *v=r1;
            return q;
        }
    }
    
    return r0;
}

// algo 1 : calcul du polynôme minimal d'une suite pour un corps Z/pZ, avec p=2, 3, 5 ou 7
int pol_min(int d, int u, int p){
    int s, t, s2, t2;
    int x=0; set(x, 1, 2*d); // x^{2d}
    int pgcd;
    int d2, deg_s, deg_t;
    int i=0, j=0, tmp;
    
    // on fait bezout tant que les degrés de t et s ne correspondent pas :
    bezout_algo1(&t, &s, u, x, p);
    
    // on faire le pgcd de t et s pour les diviser et ainsi les rendre premiers entre eux
    pgcd=bezout(&t2, &s2, t, s, p);
    
    t=div_p(t, pgcd, p);
    //s=div_p(s, pgcd, p); en commentaire, car on a pas besoin de s
    
    // on cherche à rendre rev(t) unitaire, 
    // donc le dernier chiffre différent de 0 doit être égal à 1 : 
    tmp=t;
    while((tmp%10)==0){
        tmp=tmp/10;
        i++;
    }
    tmp=tmp%10;
    
    for(j=0 ; j<p ; j++)
        if(mulmod(j, tmp, p)==1)
            break;
    
    t=mult_p(t, j, p);
    //s=mult_p(s, j, p); en commentaire, car on a pas besoin de s
    
    // on cherche le degré du polynôme minimal : 
    deg_t=degre(t);
    deg_s=degre(s);
    
    if(deg_s+1>deg_t)
        d2=deg_s+1;
    else
        d2=deg_t;
    
    // on renvoie le "reversal" de t
    return rev(t, d2);
}

int mult_mat_vect(int *A, int b, int n, int p){
    int i, j, Ai, b2, res=0;
    
    for(i=0; i<n; i++){
        Ai=A[i];
        b2=b;
        for(j=0; j<n; j++){
           set(res, addmod(get(res, n-i-1), mulmod((Ai%10), (b2%10), p), p), n-i-1);
           Ai=Ai/10;
           b2=b2/10;
        }
    }
    
    return res;
}

int mult_vect_vect(int b1, int b2, int n, int p){
    int i, j, res=0;
    
    for(i=0; i<n; i++){
        res+=mulmod(b1%10, b2%10, p);
        b1=b1/10;
        b2=b2/10;
    }
    
    return res;
}

int mult_coeff_vect(int v, int k, int p){
    int i, d=degre(v), res=0; 
    
    for(i=0; i<=d; i++){
        set(res, mulmod(v%10, k, p), i);
        v=v/10;
    }
    
    return res;
}

int add_vect_vect(int v1, int v2, int n, int p){
    int i, res=0;
    
    for(i=0; i<n; i++)
        set(res, addmod(get(v1, i), get(v2, i), p), i);
    
    return res;
}

// évaluation en a de f sur Z/pZ avec p = 2, 3, 5 ou 7
// avec la méthode de horner
int horner(int f, int a, int p){
    int i, d=degre(f), res=0;
    
    // pour avoir accès à la puissance la plus haute en premier :
    f=rev(f, d);
    
    // méthode d'Horner :
    for(i=0 ; i<=d ; i++){
        res=addmod(mulmod(res, a, p),(f%10),p);
        f=f/10;
    }
    
    return res;
}

// évaluation en a de f sur Z/pZ avec p = 2, 3, 5 ou 7
// avec la méthode de horner
int horner_mat(int f, int *A, int b, int n, int p){
    int i, d=degre(f), res=0, b_tmp;
    
    // pour avoir accès à la puissance la plus haute en premier :
    f=rev(f, d);
    
    // méthode d'Horner :
    for(i=0 ; i<=d ; i++){
        b_tmp=mult_coeff_vect(b, (f%10), p);
        res=add_p(mult_mat_vect(A, res, n, p), b_tmp, p);
        f=f/10;
    }
    
    return res;
}

// le réécrire avec horner si possible : 
int eval_pol_at_matrice(int f, int *A, int b, int n, int p){
    int i=0, d=degre(f), res=0;
    for(i=0; i<=d; i++){
        res=add_vect_vect(res, mult_coeff_vect(b, (f%10), p), n, p); // mul vect
        f=f/10;
        b=mult_mat_vect(A, b, n, p);
    }
    
    return res;
}

// algo 2 : renvoie la solution de Ay=b : 
int pol_min_matrice(int *A, int n, int b, int p){
    int i, u=100, n2, seq=0, m=1, j=0, b2=b; // mettre un u random
    if (b==0)
        return 1;
    else{
        do{
            n2=2*n;
            b2=b;
            for(i=0; i<n2; i++){
                set(seq, mult_vect_vect(u, b2, n, p), i);
                b2=mult_mat_vect(A, b2, n, p);
            }
            m=pol_min(n, seq, p);
            //cout << "m=" << m << endl;
            u=120; // mettre un u random
        }while(eval_pol_at_matrice(m, A, b, n, p)!=0 && j++<10);
        return m;
    }
}

// algo 3 : algo de Wiedemann
int Wiedemann(int *A, int n, int b, int p){
    int m = pol_min_matrice(A, n, b, p);
    int h;
    int h1, h2, x=10;
    int m0 = horner(m, 0, p);
    
    h1 = mult_p(diff_p(m, m0, p), -1+p, p);
    h2=mult_p(m0, x, p);
    
    h=div_p(h1, h2, p);
    
    return horner_mat(h, A, b, n, p);    
}

int main(int argc, char** argv) {
    //int u, v;
    
    /*cout << "bezout : " << bezout(&u, &v, 1000000, 32403, 5) << endl;
    cout << "u=" << u << " v=" << v << endl;*/
    
    cout << "pol min : " << pol_min(3, 32403, 5) << endl;
    cout << "pol min : " << pol_min(3, 10100, 3) << endl;
    
    //cout << "mult = " << mult_p(32403, 221, 5) << endl;
    //cout << "div : " << div_p(1000000, 32403, 5) << endl;
    /*cout << "div : " << div_p(1000000, 10100, 3) << endl;
    cout << "m : " << mult_p(102, 10100, 3) << endl;
    cout << "r=" << diff_p(1000000, mult_p(102, 10100, 3), 3) << endl;
    cout << "div : " << div_p(10100, 100, 3) << endl;
    cout << "m : " << mult_p(101, 100, 3) << endl;
    cout << "r=" << diff_p(10100, mult_p(100, 101, 3), 3) << endl;*/
    //cout << "div : " << div_p(1000000, 10100, 3) << endl;
    //cout << "diff = " << diff_p(1000000, 1430100, 5) << endl;
    //cout << "div : " << div_p(42, 3, 5) << endl;
    
    //cout << horner(10204, 3, 5) << endl;
    //cout << horner(2000511, 1, 7) << endl;
    //cout << horner(2000511, 2, 7) << endl;
    
    int A[3]={144, 403, 124};
    int b=312;
    //int b2=33;
    //int b3=443;
    
    //cout << "res mat fois vect : " << mult_mat_vect(A, b, 3, 5) << endl;
    //cout << "res mat fois mat : " << mult_mat_mat(A, A, 3, 5) << endl;
    //cout << "res : " << mult_mat_vect(A, b2, 3, 5) << endl;
    //cout << "res : " << mult_mat_vect(A, b3, 3, 5) << endl;
    
    cout << "Polynôme minimal : " << pol_min_matrice(A, 3, b, 5) << endl; 
    
    cout << "Wiedemann : " << Wiedemann(A, 3, b, 5) << endl;
    
    
    return 0;
}

