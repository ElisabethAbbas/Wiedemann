#include <cstdlib>
#include <iostream>
#include <tgmath.h>

using namespace std;

// Structures de données :

struct Poly {
	unsigned int* p; // contient les coefficients de p
        unsigned int corps; // corps du polynôme p
        unsigned int degre; // degré du polynôme p
        unsigned int taille; // nombre de cases
        unsigned int bits; // nombre de bits de p;
        unsigned int taille_coefficient; // nombre de bits pour un coefficient
        unsigned coefficients_par_case; // nombre de coefficients par case
};

struct Coefficient {
    unsigned int i; // numéro de ligne
    unsigned int j; // numéro de colonne
    unsigned int valeur;
};

struct Matrice{
    unsigned int lignes; // nombre de lignes
    unsigned int colonnes; // nombre de colonnes
    Coefficient *coefficients; // tableau de Coefficient
};

int get(Poly *p, unsigned int i);
int set(Poly *p, unsigned int i, unsigned int a);
Poly creer_poly_nul(int corps);
Poly creer_poly(unsigned int *t, int degre, int corps);
Poly copie_poly(Poly *p);
void print_poly (Poly *p);
void print_details_poly(Poly* p);
Poly rev(Poly p);
Poly add_poly(Poly *a, Poly *b);
Poly diff_poly(Poly *a, Poly *b);
Poly mult_poly(Poly *a, Poly *b);
Poly div_poly(Poly *a, Poly *b , Poly* r);

/*--------------------------------------------------------------------*/
/*Fonctions sur les polynômes*/
/*--------------------------------------------------------------------*/

/* crée et renvoie le polynôme p via le tableau de coefficients "t", 
 * de degré "degre" et sur le corps "corps". */
Poly creer_poly_nul(int corps){
    int b=0, j=0;
    Poly p;
    p.corps=corps;
    p.degre=0;
    
    /* On arrondit au-dessus pour trouver le nombre nécessaire
     * de bits pour un coefficient :*/
    p.taille_coefficient=ceil(log2(corps)); 
    
    // Le nombre de coefficient qu'on peut mettre par case du tableau :
    p.coefficients_par_case = (sizeof(unsigned int))*8/p.taille_coefficient;
    
    p.taille=((p.degre+1)/p.coefficients_par_case)+1;
    
    // Le nombre de bits :
    p.bits=p.taille*sizeof(unsigned int)*8;
    
    //On remplit le tableau avec des 0 :
    p.p = (unsigned int*)malloc(p.bits);
    
    for(b=0; b<p.taille; b++)
        p.p[b]=0;
        
    return p;
}

/* crée et renvoie le polynôme p via le tableau de coefficients "t", 
 * de degré "degre" et sur le corps "corps". */
Poly creer_poly(unsigned int *t, int degre, int corps){
    int i=0, j=0;
    Poly p;
    p.corps=corps;
    p.degre=degre;
    
    /* On arrondit au-dessus pour trouver le nombre nécessaire
     * de bits pour un coefficient :*/
    p.taille_coefficient=ceil(log2(corps)); 
    
    // Le nombre de coefficient qu'on peut mettre par case du tableau :
    p.coefficients_par_case = (sizeof(unsigned int))*8/p.taille_coefficient;
    
    p.taille=((p.degre+1)/p.coefficients_par_case)+1;
    
    // Le nombre de bits :
    p.bits=p.taille*sizeof(unsigned int)*8;
    
    
    //On remplit le tableau avec les coefficients du polynôme (bit-packing) :
    
    p.p = (unsigned int*)malloc(p.bits);
   
    int b=0; // numéro du bloc en cours
    
    p.p[0]=0;
    for(i=0; i<=degre; i++){
        
        // si on a remplit un bloc :
        if(j==p.coefficients_par_case){
            j=0; // on repart à 0
            b++; // on passe au bloc suivant
            p.p[b]=0;
        }
        set(&p, i, t[i]);
        j++;
    }
    
    return p;
    
}

Poly copie_poly(Poly *p){
    int b=0, j=0;
    Poly res;
    res.corps=p->corps;
    res.degre=p->degre;
    res.taille_coefficient=p->taille_coefficient; 
    res.coefficients_par_case=p->coefficients_par_case;
    res.bits=p->bits;
    res.taille=p->taille;
    
    res.p = (unsigned int*)malloc(p->bits);
    
    for(b=0; b<p->taille; b++)
        res.p[b]=p->p[b];
        
    return res;
}

// renvoie le coefficient de X^i dans le polynôme a
int get(Poly *p, unsigned int i){
    // masque permemttant de récupérer le coefficient :
    unsigned int masque=0;
    int j=0;
    
    /* le masque contient autant de 1 que de bits 
     * réservés à chaque coefficient : */
    
    do{
        masque<<=1;
        masque|=1;
    }while((++j)<p->taille_coefficient);
    
    // numéro du bloc où se trouve le coefficient cherché :
    int num_bloc=i/p->coefficients_par_case;
    
    // on récupère le block où se trouve le coefficient :
    unsigned int bloc=p->p[num_bloc];
    
    // on renvoie le coefficient : 
    return (bloc>>((i*p->taille_coefficient)%(p->coefficients_par_case*p->taille_coefficient)))&masque;
}

// modifie le coefficient de X^i dans le polynôme *p par a
int set(Poly *p, unsigned int i, unsigned int a){
    int j=0;
    int b, num_bit;
    unsigned int k;
    
    
    int num_bloc=i/p->coefficients_par_case;
    int num_coeff=i-num_bloc*p->coefficients_par_case;
    
    num_bit=num_coeff*p->taille_coefficient;
    
    /*mise à jour de la structure*/
    if((i>p->degre) && a!=0){
        p->degre=i;
        p->taille=((p->degre+1)/p->coefficients_par_case)+1;
        p->bits=(p->degre+1)*p->taille_coefficient;
        p->p=(unsigned int*)realloc(p->p, p->bits);
    }
    
    else if(i<=p->degre && (get(p,p->degre)==0)){
        k=0;
        while(p->degre>0 && get(p, p->degre)==0){
            p->degre--;
            k++;
        }
        p->taille=p->taille-k/p->taille_coefficient;
        p->bits=p->taille_coefficient*(p->degre+1);
        p->p=(unsigned int*)realloc(p->p, p->bits);
    }
    
    /* mise à jour du coefficient */
    while(j<p->taille_coefficient){
        b=(a&1);
        if (b)
            p->p[num_bloc] = p->p[num_bloc] | (1u<<(num_bit+j));
        else
            p->p[num_bloc]=p->p[num_bloc] & (~(1u<<(num_bit+j)));
        j++;
        a>>=1;
    }   
}

// affiche le polynôme *p
void print_poly (Poly *p){
    int j=0, a;
    
    /* on affiche le coefficient pour le degré 0 :*/
    // s'il existe et si le coefficient n'est pas 0
    if((j++)<=p->degre && (a=(get(p,0)))){  // on récupère le coefficient dans a
        cout << a;
    }
    
    /* on affiche le coefficient pour le degré 1 :*/
    if((j++)<=p->degre){
        if(a)
            cout << " + ";
        if(a=(get(p,1))){ // on récupère le coefficient dans a et on vérifie qu'il n'est pas nul
            if (a!=1) // s'il n'est pas égal à 1 on l'affiche
                cout << a;
            cout << "X";
        }
    }
    
    /* on affiche le reste des coefficients :*/
    while(j<=p->degre){
        if(a)
            cout << " + ";
        if(a=(get(p,j))){
            if(a!=1)
                cout << a;
            cout << "X^" << j;
        }
        
        j++;
    }    
    
    cout << endl;
}

void print_details_poly(Poly* p){
    cout << "Nombre de bits : " << p->bits << endl;
    cout << "Nombre de coefficients par case : ";
    cout << p->coefficients_par_case << endl;
    cout << "Corps : " << p->corps << endl;
    cout << "Degré : " << p->degre << endl;
    cout << "Taille : " << p->taille << endl;
    cout << "Nombre de bits d'un coefficient : "; 
    cout << p->taille_coefficient << endl;
    print_poly(p);
    cout << endl;
}

// Trouver le polynôme réciproque revP(x) = x^d P(1/x)
Poly rev(Poly p){
    unsigned int t[p.degre+1];
    for (int i = 0; i<=(p.degre); i++) {
        int a = get(&p,i);
        t[(p.degre)-i]=a;
    }
    Poly rev_p = creer_poly(t, p.degre, p.corps);
    return rev_p;
}

// Renvoie l'addition de deux polynômes à valeurs dans le même corps
Poly add_poly(Poly *a, Poly *b){
    Poly res;
    int i=0;
    
    // degré à partir duquel il faut commencer la boucle
    unsigned int d=max(a->degre, b->degre);
    
    res=creer_poly_nul(a->corps);
    for(i=d; i>=0; i--)
        set(&res, i, ((get(a, i)+get(b, i))%res.corps));
    
    return res;
}

// renvoie la différence du polynôme a par le polynôme b, 
// il faut que les deux polynômes aient des coefficients dans le même corps
Poly diff_poly(Poly *a, Poly *b){
    Poly res;
    int i=0;
    
    // degré à parti duquel il faut commencer la boucle
    unsigned int d=max(a->degre, b->degre);
    
    // tableau de coefficients nuls
    //unsigned int taille = max(a->taille, b->taille);

    res=creer_poly_nul(a->corps);
    
    //res.corps=a->corps;
    //res.degre=d;
    
    for(i=d; i>=0; i--)
        set(&res, i, ((get(a, i)+res.corps-get(b, i))%res.corps)); 
        
    /*res.coefficients_par_case=a->coefficients_par_case;
    res.taille_coefficient=a->taille_coefficient;
    
    // on met à jour res :
    unsigned int k=0;
    while(res.degre>0 && get(&res, res.degre)==0){
        res.degre--;
        k++;
    }
    res.taille=taille-k/res.taille_coefficient;
    res.bits=res.taille_coefficient*(res.degre+1);
            
    res.p=(unsigned int*)realloc(res.p, res.bits);
    */
    return res;
}

Poly mult_poly(Poly *a, Poly *b){
    Poly res;
    int r=0;
    int i,j;
    
    res=creer_poly_nul(a->corps);

    //set(&res, d, (get(a, a->degre)+1)*(get(b, (b->degre))+1)); // coeff de plus haut degré
    for(i=0; i<=a->degre; i++){
         for (j=0; j<=b->degre; j++) {
            r = get(&res, i+j); // valeur du coeff à la puissance 1O^{i+j}
            set(&res, i+j, (r+ get(a, i)*get(b, j))%res.corps);
        }
    }

    return res;
}

Poly div_poly(Poly *a, Poly *b , Poly* r){
    int kr=a->degre;
    int kb=b->degre;
    int i;
    Poly tmp, q=creer_poly_nul(a->corps);    
    
    *r=copie_poly(a);
    
    /* tant que le degré du dividende est plus grand ou égal à celui du 
     * diviseur*/
    while(kr>=kb){
        /* on cherche le coefficient qui, multiplié au diviseur, donne 
         * celui du dividende :*/
        for(i=0; i<a->corps; i++){
            if(((i*get(b, kb))%(b->corps)) == get(r, kr)){
                
                // on met à jour une variable temporaire :
                tmp=creer_poly_nul(a->corps);
                
                set(&tmp, kr-kb, i);
                // on met à jour le quotient : 
                set(&q, kr-kb, i);
                
                /* on soustrait (comme dans l'algorithme) pour obtenir le 
                 * nouveau dividende :*/
                
                tmp=mult_poly(b, &tmp);
                
                *r=diff_poly(r, &tmp);
                
                // on met à jour le degré du dividende :
                kr=r->degre;
                
                break;
            }
        }
    }    
    
    return q;
}

// comparaison de deux polynômes
int egaux (Poly* p1, Poly* p2){
    int res = 1;
    int d1 = p1->degre;
    int d2 = p2 ->degre;
    if (d1 !=d2) {
        res = 0;
    }
    else{
        for (int i = 0; i<= p1->degre; i++) {
            if (get(p1, i)!=get(p2, i)) {
                res = 0;
            }
        }
    }
    return res;
}

// ALGO DE BEZOUT
Poly bezout(Poly *u, Poly *v, Poly a, Poly b){
    int corps = u->corps;
    Poly nul = creer_poly_nul(corps);
    Poly temp =nul;
    unsigned int t0[1] = {0};
    unsigned int t1[1] = {1};
    Poly u0=creer_poly(t0, 0, corps);
    Poly u1= creer_poly(t0, 0, corps);
    Poly v0=creer_poly(t1, 0, corps);
    Poly v1=creer_poly(t1,0, corps);
    Poly ut, vt;
    Poly r0 = a;
    Poly r1 = b;
    Poly rt;
    Poly q;
    
    
    while(egaux(&r1, &nul)==0){
        q=div_poly(&r0, &r1, &nul);
        
        //cout << "q = " << q << endl;
        
        // colonne des restes :
        temp =mult_poly(&q, &r1);
        rt=diff_poly(&r0, &temp);
        r0=r1;
        r1=rt;
        temp=nul;
        
        // colonne des u :
        //cout << "u0 = " << u0 << endl;
        //cout << "u1 = " << u1 << endl;
        temp=mult_poly(&q, &u1);
        ut=diff_poly(&u0, &temp);
        u0=u1;
        u1=ut;
        temp=nul;
        
        // colonne des v :
        //cout << "v0 = " << v0 << endl;
        //cout << "v1 = " << v1 << endl;
        temp =mult_poly(&q, &v1);
        vt=diff_poly(&v0, &temp);
        v0=v1;
        v1=vt;
        temp=nul;
    }
    // mise à jour de u et v
    *u=u0;
    *v=v0;
    return r0;
}

// BEZOUT POUR ALGO 1 : tous dans le même corps
Poly bezout_algo1(Poly *u, Poly *v, Poly a, Poly b){
    int corps = u->corps;
    Poly nul = creer_poly_nul(corps);
    Poly temp =nul;
    unsigned int t0[1] = {0};
    unsigned int t1[1] = {1};
    Poly u0=creer_poly(t0, 0, corps);
    Poly u1= creer_poly(t0, 0, corps);
    Poly v0=creer_poly(t1, 0, corps);
    Poly v1=creer_poly(t1,0, corps);
    Poly ut, vt;
    Poly r0 = a;
    Poly r1 = b;
    Poly rt;
    Poly q;
    int d_n = b.degre/2;
    
    while(egaux(&r1, &nul)==0){
        q=div_poly(&r0, &r1, &nul);
        
        // colonne des restes :
        temp=mult_poly(&q, &r1);
        rt=diff_poly(&r0, &temp);
        r0=r1;
        r1=rt;
        temp=nul;
        
        // colonne des u :
        temp=mult_poly(&q, &u1);
        ut=diff_poly(&u0, &temp);
        u0=u1;
        u1=ut;
        temp=nul;
        
        // colonne des v :
        temp=mult_poly(&q, &v1);
        vt=diff_poly(&v0, &temp);
        v0=v1;
        v1=vt;
        temp=nul;
        
        // mise à jour de u et v
        if(u1.degre<=d_n && r1.degre<d_n){
            *u=u1;
            *v=r1;
            return q;
        }
    }
    return r0;
}

// ALGO 1 : calcul du polynôme minimal d'une suite
Poly pol_min(int d, Poly u){
    int corps = u.corps;
    Poly nul = creer_poly_nul(corps);
    Poly s, t, s2, t2;
    
    Poly x;
    unsigned int x2d[2*d+1];
    for (int i = 0; i<=2*d; i++) {
        x2d[i] = 0;
    }
    x2d[2*d+1] = 1;
    x = creer_poly(x2d, 2*d, corps);
    
    Poly pgcd;
    int d2, deg_s, deg_t;
    int i=0, j=0;
    
    // on fait bezout tant que les degrés de t et s ne correspondent pas :
    bezout_algo1(&t, &s, u, x);
    
    // on faire le pgcd de t et s pour les diviser et ainsi les rendre premiers entre eux
    pgcd=bezout(&t2, &s2, t, s);
    
    t=div_poly(&t, &pgcd, &nul);
    //s=div_p(s, pgcd); en commentaire, car on a pas besoin de s
    // on cherche à rendre rev(t) unitaire,
    // donc le dernier chiffre différent de 0 doit être égal à 1 :
    
    // on cherche le degré du polynôme minimal :
    deg_t=t.degre;
    deg_s=s.degre;
    
    if(deg_s+1>deg_t)
        d2=deg_s+1;
    else
        d2=deg_t;
    
    // on renvoie le "reversal" de t
    return rev(t);
}


// ------------------------------------------------------------------//
//                                 MAIN
// ------------------------------------------------------------------//

int main(int argc, char** argv) {
    Poly p, pp, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13;
    
    unsigned int tt[]={3};
    pp=creer_poly(tt, 0, 4);
    

    // TEST DE REV_P
    unsigned int t[]={2,3,4,1,0,0,1,0,4,7,5,8,6,3,6,7,12,8,3,4,4,1};
    p=creer_poly(t, 21, 23);
    cout << "P(x) = ";
    print_poly(&p);
    Poly rev_p = rev(p);
    cout << "revP(x) = ";
    print_poly(&rev_p);
    cout << endl;
    
    /*
    // TEST ALGO 1:
    unsigned int u[] = {3,0,4,2,3,0};
    U = creer_poly(u, 5, 5);
    p = pol_min(3, U);
    print_poly(&p);
    */
    
    /*
    // TESTS DE DIFF_P
    unsigned int t2[]={2,4,3};
    p2=creer_poly(t2,2,7);
    unsigned int t3[]={1,4,6};
    p3=creer_poly(t3,2,7);
    p4=diff_poly(&p2, &p3);
    cout << "DIFFERENCE DANS Z/7Z :" << endl;
    print_poly(&p2);
    cout << "-" << endl;
    print_poly(&p3);
    cout << "=" << endl;
    print_poly(&p4);
    cout << endl;
    
    unsigned int t5[]={1,1,4,2,0,1};
    unsigned int t6[]={1,3,4,1,0,2,1,1,3};
    p5=creer_poly(t5, 5, 5);
    p6=creer_poly(t6, 8, 5);
    p7=diff_poly(&p5, &p6);
    cout << "DIFFERENCE DANS Z/5Z :" << endl;
    print_poly(&p5);
    cout << "-" << endl;
    print_poly(&p6);
    cout << "= " << endl;
    print_poly(&p7);
    cout << endl;
    
    
    // TESTS DE ADD_P
    
    cout << "ADDITION DANS Z/7Z :" << endl;
    p8=add_poly(&p2, &p3);
    print_poly(&p2);
    cout << "-" << endl;
    print_poly(&p3);
    cout << "=" << endl;
    print_poly(&p8);
    cout << endl;
    
    cout << "ADDITION DANS Z/5Z : " << endl;
    p9 = add_poly(&p5, &p6);
    print_poly(&p5);
    cout << "+ " << endl;
    print_poly(&p6);
    cout << "= " << endl;
    print_poly(&p9);
    cout << endl;
    
    // TEST MULTIPLICATION
    cout << "MULTIPLICATION DANS Z/7Z : "<<endl;
    p10 = mult_poly(&p2, &p3);
    print_poly(&p2);
    cout << "x" << endl;
    print_poly(&p3);
    cout << "=" << endl;
    print_poly(&p10);
    cout << endl;
    
    
    unsigned int t11[]={1,4,1};
    unsigned int t12[]={1,1};
    p11=creer_poly(t11, 2, 11);
    p12=creer_poly(t12, 1, 11);
    p13=mult_poly(&p11, &p12);
    cout << "MULTIPLICATION DANS Z/11Z : "<<endl;
    print_poly(&p11);
    cout << "x" << endl;
    print_poly(&p12);
    cout << "=" << endl;
    print_poly(&p13);
    cout << endl;
    */
    
    return 0;
}

