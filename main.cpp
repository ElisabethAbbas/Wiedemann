#include <cstdlib>
#include <iostream>
#include <tgmath.h>

#define TAUXREMPLISSAGE 10
#define TAILLE 200
#define CORPS 5


using namespace std;


struct Poly {
	unsigned int* p; // contient les coefficients de p
        unsigned int corps; // corps du polynôme p
        unsigned int degre; // degré du polynôme p
        unsigned int taille; // nombre de cases
        unsigned int bits; // nombre de bits de p;
        unsigned int taille_coefficient; // nombre de bits pour un coefficient
        unsigned coefficients_par_case; // nombre de coefficients par case
};

struct Vecteur {
    unsigned int n; // taille du vecteur
    unsigned int taille; // nombre de cases
    unsigned int nb_coefficients; // nombre de coefficients
    unsigned int *coefficients; // tableau de coefficients non nuls
    unsigned int corps;
    unsigned int taille_coefficient; // nombre de bits d'un coefficient
    unsigned int taille_n; // nombre de bits pour positionner un coefficient
    unsigned int coefficients_par_case; // nombre de coefficients par case
};

struct Matrice {
    unsigned int n; // taille du vecteur
    unsigned int taille; // nombre de cases
    unsigned int nb_coefficients; // nombre de coefficients
    unsigned int *coefficients; // tableau de coefficients non nuls
    unsigned int corps;
    unsigned int taille_coefficient; // nombre de bits d'un coefficient
    unsigned int taille_n; // nombre de bits nécessaires pour i ou j
    unsigned int coefficients_par_case; // nombre de coefficients par case
};

unsigned int addmod(unsigned int a, unsigned int b, unsigned int p);
unsigned int mulmod(unsigned int a, unsigned int b, unsigned int p);

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

int egaux (Poly* p1, Poly* p2);

Vecteur creer_vecteur(unsigned int t[][2], unsigned int taille, 
        unsigned int n, unsigned int corps);
void set(Vecteur *v, unsigned int i, unsigned int a);
unsigned int get(Vecteur *v, unsigned int i);
void print_details_vecteur(Vecteur *v);
void print_vecteur(Vecteur *v);

Matrice creer_matrice(unsigned int t[][3], unsigned int taille, 
        unsigned int n, unsigned int corps);
void set(Matrice *m, unsigned int i, unsigned int j, unsigned int a);
unsigned int get(Matrice *m, unsigned int i, unsigned int j);
void print_matrice(Matrice *m);


Vecteur add_vect_vect(Vecteur *v, Vecteur *v2);
Vecteur mult_coeff_vect(Vecteur *v, int k);
unsigned int mult_vect_vect(Vecteur *v, Vecteur *b2);

Vecteur horner_mat(Poly f, Matrice *A, Vecteur b);

Poly vecteur_to_poly(Vecteur *v);


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
    
    for(b=0; b<p->taille; b++){
        res.p[b]=p->p[b];
    }
        
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
    return (bloc>>((i*p->taille_coefficient)%(p->coefficients_par_case
            *p->taille_coefficient)))&masque;
}

// modifie le coefficient de X^i dans le polynôme *p par a
int set(Poly *p, unsigned int i, unsigned int a){
    int j=0, l;
    int b, num_bit;
    unsigned int k;
    unsigned int tmp;
    
    
    int num_bloc=i/p->coefficients_par_case;
    int num_coeff=i-num_bloc*p->coefficients_par_case;
    
    num_bit=num_coeff*p->taille_coefficient;
    
    
    /*mise à jour de la structure*/
    if((i>p->degre) && a!=0){
        p->degre=i;
        tmp=p->taille;
        p->taille=((p->degre+1)/p->coefficients_par_case)+1; 
        for(l=tmp; l<p->taille; l++)
            p->p[l]=0;
        p->bits=(p->degre+1)*p->taille_coefficient;
        p->p=(unsigned int*)realloc(p->p, p->bits);
    }
    
    else if(i<=p->degre && (get(p,p->degre)==0) && a==0){ //
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
    if(((j++)<=p->degre && (a=(get(p,0)))) || p->degre==0){  
        // on récupère le coefficient dans a
        cout << a;
    }
    
    /* on affiche le coefficient pour le degré 1 :*/
    if((j++)<=p->degre){
        if(a)
            cout << " + ";
        if(a=(get(p,1))){ // on récupère le coefficient dans a et 
            // on vérifie qu'il n'est pas nul
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

    res=creer_poly_nul(a->corps);
    
    for(i=d; i>=0; i--)
        set(&res, i, ((get(a, i)+res.corps-get(b, i))%res.corps)); 
        
    
    return res;
}

Poly mult_poly(Poly *a, Poly *b){
    Poly res;
    int r=0;
    int i,j;
    
    res=creer_poly_nul(a->corps);

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
    
    Poly p_nul=creer_poly_nul(a->corps);
    
    *r=copie_poly(a);
    
    /* tant que le degré du dividende est plus grand ou égal à celui du 
     * diviseur*/
    while(kr>=kb && egaux(r, &p_nul)!=1){ // r n'est pas nul !
        
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

// ------------------------------------------------------------------//
//                             Vecteurs
// ------------------------------------------------------------------//

Vecteur creer_vecteur(unsigned int t[][2], unsigned int taille, 
        unsigned int n, unsigned int corps){
    Vecteur v;
    unsigned int i, b=0;
    
    v.n=n;
    v.nb_coefficients=1;
    v.corps=corps;
    v.taille_coefficient=ceil(log2(v.corps));
    v.taille_n=ceil(log2(v.n));
    v.coefficients_par_case=sizeof(unsigned int)*
            8/(v.taille_coefficient+v.taille_n);
    //v.taille=(v.nb_coefficients/v.coefficients_par_case)+1;
    v.taille=1;
    v.coefficients=(unsigned int*)malloc(v.nb_coefficients*
            (v.taille_coefficient+v.taille_n));
    
    set(&v, 0, 0);
    
    for(i=0; i<taille; i++){
        if(i%v.coefficients_par_case==0){
            v.coefficients[b]=0;
            b++;
        }
        set(&v, t[i][0], t[i][1]);
    }
    return v;
}

Vecteur creer_vecteur_nul(unsigned int n, unsigned int corps){
    Vecteur v;
    
    v.n=n;
    v.nb_coefficients=1;
    v.corps=corps;
    v.taille_coefficient=ceil(log2(v.corps));
    v.taille_n=ceil(log2(v.n));
    v.coefficients_par_case=sizeof(unsigned int)*8/(v.taille_coefficient
            +v.taille_n);
    v.taille=1;
    
    v.coefficients=(unsigned int*)malloc(v.nb_coefficients*
            (v.taille_coefficient+v.taille_n));
    
    v.coefficients[0]=0;
    
    set(&v, 0, 0);
    
    return v;
}

void set(Vecteur *v, unsigned int i, unsigned int a){
    unsigned int j=0, b, k;
    int num_bit, bloc=0, num_bloc;
    unsigned int masque_c=0, masque_n=0;
    
    j=0;
    do{
        masque_n<<=1;
        masque_n|=1;
    }while((++j)<v->taille_n);
    
    j=0;
    do{
        masque_c<<=1;
        masque_c|=1;
    }while((++j)<v->taille_coefficient);
    
    // on cherche si le coefficient existe déjà
    for(num_bloc=0; num_bloc<v->taille; num_bloc++){
        bloc=v->coefficients[num_bloc];
        for(k=0; k<v->coefficients_par_case; k++){
            if((bloc&masque_n)==i){
                num_bit=(k*(v->taille_n+v->taille_coefficient));
                num_bit+=v->taille_n;
                j=0;
                while(j<v->taille_coefficient){
                    b=(a&1);
                    if (b)
                        v->coefficients[num_bloc] = v->coefficients[num_bloc] 
                                | (1u<<(num_bit+j));
                    else
                        v->coefficients[num_bloc] = v->coefficients[num_bloc] 
                                & (~(1u<<(num_bit+j)));
                    j++;
                    a>>=1;
                }
                return;
            }
            bloc>>=(v->taille_n+v->taille_coefficient);
        }
    }
    
    
    
    if((v->taille*v->coefficients_par_case)==(v->nb_coefficients)){
        num_bit=0;
        num_bloc=v->taille;
        v->taille++;
        v->nb_coefficients++;
        v->coefficients=(unsigned int*)realloc(v->coefficients, 
                v->taille*(v->coefficients_par_case*
                (v->taille_coefficient+v->taille_n)));
    }
    else{
        num_bit=(v->nb_coefficients%v->coefficients_par_case)*
                (v->taille_coefficient+v->taille_n);
        num_bloc=v->taille-1;
        v->nb_coefficients++;
        v->coefficients=(unsigned int*)realloc(v->coefficients, v->taille*
                (v->coefficients_par_case*(v->taille_coefficient
                +v->taille_n)));
    }
    
    j=0;
    while(j<v->taille_n){
        b=(i&1);
        if (b)
            v->coefficients[num_bloc] = v->coefficients[num_bloc] | 
                    (1u<<(num_bit+j));
        else
            v->coefficients[num_bloc] = v->coefficients[num_bloc] & 
                    (~(1u<<(num_bit+j)));
        j++;
        i>>=1;
    }
    
    j=0;
    num_bit+=v->taille_n;
    
    while(j<v->taille_coefficient){
        b=(a&1);
        if (b)
            v->coefficients[num_bloc] = v->coefficients[num_bloc] | 
                    (1u<<(num_bit+j));
        else
            v->coefficients[num_bloc] = v->coefficients[num_bloc] & 
                    (~(1u<<(num_bit+j)));
        j++;
        a>>=1;
    }
    
    return;
}

Vecteur copie_vecteur(Vecteur *p){
    int b=0, j=0;
    Vecteur res;
    res.corps=p->corps;
    res.n=p->n;
    res.nb_coefficients=p->nb_coefficients;
    res.taille_n=p->taille_n;
    res.taille_coefficient=p->taille_coefficient; 
    res.coefficients_par_case=p->coefficients_par_case;
    res.taille=p->taille;
    
    res.coefficients = (unsigned int*)malloc((p->taille_coefficient+p->taille_n)*p->nb_coefficients);
    
    for(b=0; b<p->taille; b++)
        res.coefficients[b]=p->coefficients[b];
        
    return res;
}

Matrice creer_matrice(unsigned int t[][3], unsigned int taille, 
        unsigned int n, unsigned int corps){
    Matrice m;
    unsigned int i, b=0;
    
    m.n=n;
    m.nb_coefficients=1;
    m.corps=corps;
    m.taille_coefficient=ceil(log2(m.corps));
    m.taille_n=ceil(log2(m.n));
    m.coefficients_par_case=sizeof(unsigned int)*8/(m.taille_coefficient+
            (m.taille_n*2));
    m.taille=1;
    m.coefficients=(unsigned int*)malloc(m.nb_coefficients*
            (m.taille_coefficient+(m.taille_n*2)));
    
    m.coefficients[0]=0;
    set(&m, 0, 0, 0);
    
    for(i=0; i<taille; i++){
        if(i%m.coefficients_par_case==0){
            m.coefficients[b]=0;
            b++;
        }
        set(&m, t[i][0], t[i][1], t[i][2]);
    }
    
    return m;
}

Matrice creer_matrice_nulle(unsigned int n, unsigned int corps){
    Matrice m;
    unsigned int i, b=0;
    
    m.n=n;
    m.nb_coefficients=1;
    m.corps=corps;
    m.taille_coefficient=ceil(log2(m.corps));
    m.taille_n=ceil(log2(m.n));
    m.coefficients_par_case=sizeof(unsigned int)*8/(m.taille_coefficient+
            (m.taille_n*2));
    m.taille=1;
    m.coefficients=(unsigned int*)malloc(m.nb_coefficients*
            (m.taille_coefficient+(m.taille_n*2)));
    
    m.coefficients[0]=0;
    set(&m, 0, 0, 0);
    
    return m;
}

// pour la création de grandes matrices
void set2(Matrice *m, unsigned int i, unsigned int j, unsigned int a){
    unsigned int b, k=0, k2=0;
    int num_bit, bloc=0, num_bloc;
    unsigned int masque_c=0, masque_n=0;
    
    // s'il faut rajouter une case au tableau :
    if((m->taille*m->coefficients_par_case)==(m->nb_coefficients)){
        num_bit=0;
        num_bloc=m->taille;
        m->taille++;
        m->nb_coefficients++;
        m->coefficients=(unsigned int*)realloc(m->coefficients, 
                m->taille*(m->coefficients_par_case*
                (m->taille_coefficient+(m->taille_n*2))));
    }
    // sinon :
    else{
        num_bit=(m->nb_coefficients%m->coefficients_par_case)*
                (m->taille_coefficient+(m->taille_n*2));
        num_bloc=m->taille-1;
        m->nb_coefficients++;
        m->coefficients=(unsigned int*)realloc(m->coefficients, 
                m->taille*(m->coefficients_par_case*
                (m->taille_coefficient+(m->taille_n*2))));
    }
    
    // On ajoute la ligne i :
    k2=0;
    while(k2<m->taille_n){
        b=(i&1);
        if (b)
            m->coefficients[num_bloc] = m->coefficients[num_bloc] | 
                    (1u<<(num_bit+k2));
        else
            m->coefficients[num_bloc] = m->coefficients[num_bloc] & 
                    (~(1u<<(num_bit+k2)));
        k2++;
        i>>=1;
    }
    
    num_bit+=m->taille_n;
    
    // On ajoute la colonne j
    k2=0;
    while(k2<m->taille_n){
        b=(j&1);
        if (b)
            m->coefficients[num_bloc] = m->coefficients[num_bloc] | 
                    (1u<<(num_bit+k2));
        else
            m->coefficients[num_bloc] = m->coefficients[num_bloc] & 
                    (~(1u<<(num_bit+k2)));
        k2++;
        j>>=1;
    }
    
    num_bit+=m->taille_n;
    
    // On ajoute le coefficient
    k2=0; 
    while(k2<m->taille_coefficient){
        b=(a&1);
        if (b)
            m->coefficients[num_bloc] = m->coefficients[num_bloc] | 
                    (1u<<(num_bit+k2));
        else
            m->coefficients[num_bloc] = m->coefficients[num_bloc] & 
                    (~(1u<<(num_bit+k2)));
        k2++;
        a>>=1;
    }
}

void set(Matrice *m, unsigned int i, unsigned int j, unsigned int a){
    unsigned int b, k=0, k2=0;
    int num_bit, bloc=0, num_bloc;
    unsigned int masque_c=0, masque_n=0;
    
    k2=0;
    do{
        masque_n<<=1;
        masque_n|=1;
    }while((++k2)<m->taille_n);
    
    k2=0;
    do{
        masque_c<<=1;
        masque_c|=1;
    }while((++k2)<m->taille_coefficient);
    
    // on cherche si le coefficient existe déjà
    for(num_bloc=0; num_bloc<m->taille; num_bloc++){
        bloc=m->coefficients[num_bloc];
        for(k=0; k<m->coefficients_par_case; k++){
            if((bloc&masque_n)==i && ((bloc>>m->taille_n)&masque_n)==j){
                num_bit=(k*((m->taille_n*2)+m->taille_coefficient));
                num_bit+=(m->taille_n*2);
                k2=0;
                while(k2<m->taille_coefficient){
                    b=(a&1);
                    if (b)
                        m->coefficients[num_bloc] = 
                                m->coefficients[num_bloc] | 
                                (1u<<(num_bit+k2));
                    else
                        m->coefficients[num_bloc] = 
                                m->coefficients[num_bloc] & 
                                (~(1u<<(num_bit+k2)));
                    k2++;
                    a>>=1;
                }
                return;
            }
            bloc>>=((m->taille_n*2)+m->taille_coefficient);
        }
    }
    
    // s'il faut rajouter une case au tableau :
    if((m->taille*m->coefficients_par_case)==(m->nb_coefficients)){
        num_bit=0;
        num_bloc=m->taille;
        m->taille++;
        m->nb_coefficients++;
        m->coefficients=(unsigned int*)realloc(m->coefficients, 
                m->taille*(m->coefficients_par_case*
                (m->taille_coefficient+(m->taille_n*2))));
    }
    // sinon :
    else{
        num_bit=(m->nb_coefficients%m->coefficients_par_case)*
                (m->taille_coefficient+(m->taille_n*2));
        num_bloc=m->taille-1;
        m->nb_coefficients++;
        m->coefficients=(unsigned int*)realloc(m->coefficients, 
                m->taille*(m->coefficients_par_case*
                (m->taille_coefficient+(m->taille_n*2))));
    }
    
    // On ajoute la ligne i :
    k2=0;
    while(k2<m->taille_n){
        b=(i&1);
        if (b)
            m->coefficients[num_bloc] = m->coefficients[num_bloc] | 
                    (1u<<(num_bit+k2));
        else
            m->coefficients[num_bloc] = m->coefficients[num_bloc] & 
                    (~(1u<<(num_bit+k2)));
        k2++;
        i>>=1;
    }
    
    num_bit+=m->taille_n;
    
    // On ajoute la colonne j
    k2=0;
    while(k2<m->taille_n){
        b=(j&1);
        if (b)
            m->coefficients[num_bloc] = m->coefficients[num_bloc] | 
                    (1u<<(num_bit+k2));
        else
            m->coefficients[num_bloc] = m->coefficients[num_bloc] & 
                    (~(1u<<(num_bit+k2)));
        k2++;
        j>>=1;
    }
    
    num_bit+=m->taille_n;
    
    // On ajoute le coefficient
    k2=0; 
    while(k2<m->taille_coefficient){
        b=(a&1);
        if (b)
            m->coefficients[num_bloc] = m->coefficients[num_bloc] | 
                    (1u<<(num_bit+k2));
        else
            m->coefficients[num_bloc] = m->coefficients[num_bloc] & 
                    (~(1u<<(num_bit+k2)));
        k2++;
        a>>=1;
    }
    
    return;
}

unsigned int get(Vecteur *v, unsigned int i){
    // masque permemttant de récupérer le coefficient :
    unsigned int masque_n=0, masque_c=0;
    int j,b;
    unsigned int bloc;
    
    /* le masque contient autant de 1 que de bits 
     * réservés à chaque coefficient : */
    
    j=0;
    do{
        masque_n<<=1;
        masque_n|=1;
    }while((++j)<v->taille_n);
    
    j=0;
    do{
        masque_c<<=1;
        masque_c|=1;
    }while((++j)<v->taille_coefficient);
    
    for(b=0; b<v->taille; b++){
        bloc=v->coefficients[b];
        for(j=0; j<v->coefficients_par_case; j++){
            if((bloc&masque_n)==i)
                return (bloc>>v->taille_n)&masque_c;
            bloc>>=(v->taille_n+v->taille_coefficient);
        }
    }
    
    return 0;
}

unsigned int get(Matrice *m, unsigned int i, unsigned int j){
    // masque permemttant de récupérer le coefficient :
    unsigned int masque_n=0, masque_c=0;
    int k,b;
    unsigned int bloc;
    
    /* le masque contient autant de 1 que de bits 
     * réservés à chaque coefficient : */
    
    k=0;
    do{
        masque_n<<=1;
        masque_n|=1;
    }while((++k)<m->taille_n);
    
    k=0;
    do{
        masque_c<<=1;
        masque_c|=1;
    }while((++k)<m->taille_coefficient);
    
    for(b=0; b<m->taille; b++){
        bloc=m->coefficients[b];
        for(k=0; k<m->coefficients_par_case; k++){
            if((bloc&masque_n)==i && ((bloc>>m->taille_n)&masque_n)==j)
                return (bloc>>(m->taille_n*2))&masque_c;
            bloc>>=((m->taille_n*2)+m->taille_coefficient);
        }
    }
    
    return 0;
}


void print_vecteur(Vecteur *v){
    int i;
    cout << "(";
    for(i=0; i<v->n; i++)
        cout << get(v, i) << " ";
    cout << ")";
    cout << endl;
}

void print_matrice(Matrice *m){
    int i,j;
    
    for(i=0; i<m->n; i++){
        cout << " | ";
        for(j=0; j<m->n; j++)
            cout << get(m, i, j) << " ";
        cout << " | " << endl;
    }
    
    cout << endl;
}

// Multiplication de deux vecteurs (même taille et même corps) :
unsigned int mult_vect_vect(Vecteur *v1, Vecteur *v2){
    int i;
    unsigned int res=0;
    
    Vecteur *v, *v_2;
    
    if(v1->nb_coefficients>v2->nb_coefficients){
        v=v1; v_2=v2;
    }
    else{
        v=v2; v_2=v1; 
    }
    
    // masque permemttant de récupérer le coefficient :
    unsigned int masque_n=0, masque_c=0;
    int j,b;
    unsigned int bloc;
    
    /* le masque contient autant de 1 que de bits 
     * réservés à chaque coefficient : */
    
    j=0;
    do{
        masque_n<<=1;
        masque_n|=1;
    }while((++j)<v->taille_n);
    
    j=0;
    do{
        masque_c<<=1;
        masque_c|=1;
    }while((++j)<v->taille_coefficient);
    
    for(b=0; b<v->taille; b++){
        bloc=v->coefficients[b];
        for(j=0; j<v->coefficients_par_case; j++){
            res=addmod(res, mulmod(((bloc>>v->taille_n)&masque_c), get(v_2, 
                    (bloc&masque_n)), v->corps), v->corps);
            bloc>>=(v->taille_n+v->taille_coefficient);
        }
    }
    
    return res%v->corps;
    
    /*for(i=0; i<b1->n; i++){
        res+=(get(b1, i)*get(b2, i))%b1->corps;
    }
    
    return res%b1->corps; */
}

Vecteur mult_coeff_vect(Vecteur *v, int k){
    Vecteur res = creer_vecteur_nul(v->n, v->corps);
    unsigned int i, c;
    for(i=0; i<v->n; i++){
        c=get(v, i);
        if(c!=0)
            set(&res, i, mulmod(c,k, v->corps));
    }
    
    return res;
}

Vecteur add_vect_vect(Vecteur *v1, Vecteur *v2){
    //Vecteur *v, *v_2;
    int i;
    Vecteur res=creer_vecteur_nul(v1->n, v1->corps);
    
    for(i=0; i<v1->n; i++)
        set(&res, i, addmod(get(v1, i), get(v2, i), res.corps));
    
    return res;
    /*if(v1->nb_coefficients>v2->nb_coefficients){
        v=v1; v_2=v2;
    }
    else{
        v=v2; v_2=v1;
    }
    
    // masque permemttant de récupérer le coefficient :
    unsigned int masque_n=0, masque_c=0;
    int j,b;
    unsigned int bloc, itmp;*/
    
    /* le masque contient autant de 1 que de bits 
     * réservés à chaque coefficient : */
    
    /*j=0;
    do{
        masque_n<<=1;
        masque_n|=1;
    }while((++j)<v->taille_n);
    
    j=0;
    do{
        masque_c<<=1;
        masque_c|=1;
    }while((++j)<v->taille_coefficient);
    
    for(b=0; b<v->taille; b++){
        bloc=v->coefficients[b];
        for(j=0; j<v->coefficients_par_case; j++){
            itmp=(bloc&masque_n);
            set(&res, itmp, addmod(((bloc>>v->taille_n)&masque_c), 
                    get(v_2, itmp), v->corps));
            bloc>>=(v->taille_n+v->taille_coefficient);
        }
    }
    
    return res;*/
    
    
    
    /*for(i=0; i<n; i++)
        set(res, addmod(get(v1, i), get(v2, i), p), i);
    
    return res;*/
}

Vecteur mult_mat_vect(Matrice *m, Vecteur *v){
    Vecteur res=creer_vecteur_nul(v->n, m->corps);
    
    
    // masque permemttant de récupérer le coefficient :
    unsigned int masque_n=0, masque_c=0;
    unsigned int k,b,l,c,r1;
    unsigned int bloc, jtmp;
    
    /* le masque contient autant de 1 que de bits 
     * réservés à chaque coefficient : */
    
    k=0;
    do{
        masque_n<<=1;
        masque_n|=1;
    }while((++k)<m->taille_n);
    
    k=0;
    do{
        masque_c<<=1;
        masque_c|=1;
    }while((++k)<m->taille_coefficient);
    
    
    
    for(l=0; l<m->n; l++){
        //c=0;
        for(b=0; b<v->taille; b++){
            bloc=v->coefficients[b];
            for(k=0; k<v->coefficients_par_case; k++){
                jtmp=(bloc&masque_n);
                r1= mulmod(get(m, l, jtmp), ((bloc>>v->taille_n)&masque_c), 
                        m->corps);
                if(r1!=0)
                    set(&res, l, addmod(get(&res, l), r1, v->corps));
                bloc>>=((v->taille_n)+v->taille_coefficient);
            }
        }
    }
    
    return res;
}

/*--------------------------------------------------------------------*/
/*Fonctions sur des nombres modulo p */
/*--------------------------------------------------------------------*/

// (a+b) mod p
unsigned int addmod(unsigned int a, unsigned int b, unsigned int p){
    return (a+b)%p;
}

// (a*b) mod p
unsigned int mulmod(unsigned int a, unsigned int b, unsigned int p){
    return (a*b)%p;
}

void print_details_vecteur(Vecteur *v){
    cout << "Nombre de coefficients par case : ";
    cout << v->coefficients_par_case << endl;
    cout << "Corps : " << v->corps << endl;
    cout << "Nb de coefficients : " << v->nb_coefficients << endl;
    cout << "Nb de cases : " << v->taille << endl;
    cout << "Nombre de bits d'un coefficient : "; 
    cout << v->taille_coefficient << endl;
    cout << "Nombre de bits d'un indice : "; 
    cout << v->taille_n << endl;
    print_vecteur(v);
    cout << endl;
}

int egaux (Poly* p1, Poly* p2){
    int d1 = p1->degre;
    int d2 = p2->degre;
    if (d1 !=d2) {
        return 0;
    }
    else{
        for (int i = 0; i<= p1->degre; i++) {
            if (get(p1, i)!=get(p2, i)) {
                return 0;
            }
        }
    }
    
    return 1;
}

int egaux (Vecteur* p1, Vecteur* p2){
    int d1 = p1->n;
    int d2 = p2->n;
    int dmax;
    int i;
    if (d1 !=d2) {
        return 0;
    }
    else{
        dmax=max(p1->taille, p2->taille);
        for (i = 0; i <= dmax; i++) {
            if (get(p1, i)!=get(p2, i)) {
                return 0;
            }
        }
    }
    
    return 1;
}

// renvoie le pgcd et met à jour les coefficients u et v, 
// sur les polynômes a et b et dans un corps Z/pZ, avec p=2, 3, 5 ou 7
Poly bezout(Poly *u, Poly *v, Poly a, Poly b){
    
    Poly u0=creer_poly_nul(a.corps);
    set(&u0, 0, 1);
    Poly v0=creer_poly_nul(a.corps);
    Poly u1=creer_poly_nul(a.corps);
    Poly v1=creer_poly_nul(a.corps);
    set(&v1, 0, 1);
    
    Poly ut=creer_poly_nul(a.corps);
    Poly vt=creer_poly_nul(a.corps);
    
    Poly r0=copie_poly(&a);
    Poly r1=copie_poly(&b);
    Poly rt=creer_poly_nul(a.corps);
    Poly q=creer_poly_nul(a.corps);
    Poly p_nul=creer_poly_nul(a.corps);
    
    Poly tmp=creer_poly_nul(a.corps);
    
    while(egaux(&r1,&p_nul)!=1){
        q=div_poly(&r0, &r1, &tmp);
        
        //cout << "q = " << q << endl;
        
        // colonne des restes :
        tmp=mult_poly(&q, &r1);
        rt=diff_poly(&r0, &tmp);
        r0=copie_poly(&r1); 
        r1=copie_poly(&rt);
        
        // colonne des u : 
        //cout << "u0 = " << u0 << endl;
        //cout << "u1 = " << u1 << endl;
        tmp=mult_poly(&q, &u1);
        ut=diff_poly(&u0, &tmp);
        u0=copie_poly(&u1);
        u1=copie_poly(&ut);

        // colonne des v :
        //cout << "v0 = " << v0 << endl;
        //cout << "v1 = " << v1 << endl;
        tmp=mult_poly(&q, &v1);
        vt=diff_poly(&v0, &tmp);
        v0=copie_poly(&v1);
        v1=copie_poly(&vt);
    }
    
    // mise à jour de u et v
    *u=u0;    
    *v=v0;
    
    return r0;
}

// renvoie le pgcd et met à jour les coefficients u et v, 
// sur les polynômes a et b et dans un corps Z/pZ, avec p=2, 3, 5 ou 7
Poly bezout_algo1(Poly *u, Poly *v, Poly a, Poly b){
    Poly u0=creer_poly_nul(a.corps);
    set(&u0, 0, 1);
    Poly v0=creer_poly_nul(a.corps);
    Poly u1=creer_poly_nul(a.corps);
    Poly v1=creer_poly_nul(a.corps);
    set(&v1, 0, 1);
    
    Poly ut=creer_poly_nul(a.corps);
    Poly vt=creer_poly_nul(a.corps);
    
    Poly r0=copie_poly(&a);
    Poly r1=copie_poly(&b);
    
    Poly rt=creer_poly_nul(a.corps);
    Poly q=creer_poly_nul(a.corps);
    Poly p_nul=creer_poly_nul(a.corps);
    
    Poly tmp=creer_poly_nul(a.corps);
    
    
    int d_n = b.degre/2;
    
    
    while(egaux(&r1,&p_nul)!=1){
        q=div_poly(&r0, &r1, &tmp);
        //cout << "q = " << q << endl;
        
        // colonne des restes :
        tmp=mult_poly(&q, &r1);
        rt=diff_poly(&r0, &tmp);
        /*r0=copie_poly(&r1); 
        r1=copie_poly(&rt);*/
        r0=r1; 
        r1=rt;
        
        // colonne des u : 
        //cout << "u0 = " << u0 << endl;
        //cout << "u1 = " << u1 << endl;
        tmp=mult_poly(&q, &u1);
        ut=diff_poly(&u0, &tmp);
        /*u0=copie_poly(&u1);
        u1=copie_poly(&ut);*/
        u0=u1;
        u1=ut;
        
        // colonne des v :
        //cout << "v0 = " << v0 << endl;
        //cout << "v1 = " << v1 << endl;
        tmp=mult_poly(&q, &v1);
        vt=diff_poly(&v0, &tmp);
        /*v0=copie_poly(&v1);
        v1=copie_poly(&vt);*/
        v0=v1;
        v1=vt;
        
        if(u1.degre<=d_n && r1.degre<d_n){
            *u=u1;    
            *v=r1;
            return q;
        }
    }
    
    // mise à jour de u et v
    *u=u0;    
    *v=v0;
    
    return r0;
}

// évaluation en a de f sur Z/pZ avec p = 2, 3, 5 ou 7
// avec la méthode de horner
Vecteur horner_mat(Poly f, Matrice *A, Vecteur b){
    int i;
    unsigned int j, d=f.degre, c;
    Vecteur res=creer_vecteur_nul(b.n, b.corps);
    Vecteur b_tmp, b_tmp2;
    
    // méthode d'Horner :
    for(i=d ; i>=0 ; i--){
        b_tmp=mult_coeff_vect(&b, get(&f, i));
        b_tmp2=mult_mat_vect(A, &res);
        res=add_vect_vect(&b_tmp2, &b_tmp);
    }
    
    return res;
}


// renvoie le pgcd et met à jour les coefficients u et v, 
// sur les polynômes a et b et dans un corps Z/pZ, avec p=2, 3, 5 ou 7
// avec les conditions de l'algo 1
/*int bezout_algo1(int *u, int *v, int a, int b, int p){
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
}*/

// ALGO 1 : calcul du polynôme minimal d'une suite
Poly pol_min(int d, Poly u){
    int corps = u.corps;
    Poly nul = creer_poly_nul(corps);
    Poly s, t, s2, t2;
    Poly x;
    Poly pgcd;
    int d2, deg_s, deg_t;
    int i=0, j=0;
    
    unsigned int x2d[2*d+1];
    for (int i = 0; i<2*d; i++) {
        x2d[i] = 0;
    }
    
    x2d[2*d] = 1;    
    
    x = creer_poly(x2d, 2*d, corps);    
    for(i=0; i<x.taille-1; i++)
        x.p[i]=0;
    
    
    // on fait bezout tant que les degrés de t et s ne correspondent pas :
    bezout_algo1(&t, &s, u, x);
    
    // on fait le pgcd de t et s pour les diviser et ainsi les rendre premiers entre eux
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

// algo 2 : renvoie le polmin de (A^ib): 
Poly pol_min_matrice(Matrice *A, Vecteur b){
    int i, n2, j=0, r, tmp2;
    Vecteur b2 = copie_vecteur(&b);
    int compteur = 0;
    Vecteur u = creer_vecteur_nul(b.n, b.corps);
    Vecteur v_nul=creer_vecteur_nul(b.n, b.corps);
    
    Poly m=creer_poly_nul(b.corps);
    set(&m, 0, 1);
    
    Vecteur seq=creer_vecteur_nul(b.n, b.corps);
    Vecteur tmp;
    
    Poly p_seq;
    
    if (egaux(&b,&v_nul)==1)
        return m;
    else{
        tmp=horner_mat(m, A, b);
        while(egaux(&tmp, &v_nul)!=1 &&(j++ <=10)) {
            u=creer_vecteur_nul(b.n, b.corps);
            for (i = 0; i<b.n; i++) {
                r = rand()%b.corps;
                set(&u, i, r);
            }
            compteur++;
            
            cout << compteur << ". u = ";
            print_vecteur(&u);
                        
            n2=2*b.n;
            b2=copie_vecteur(&b);
            
            seq=creer_vecteur_nul(2*b.n, b.corps);
            for(i=0; i<n2; i++){
                set(&seq, i, mult_vect_vect(&u, &b2));
                b2=mult_mat_vect(A, &b2);
            }
            
            p_seq=vecteur_to_poly(&seq);
            
            m=pol_min(b.n, p_seq);
            
            tmp=horner_mat(m, A, b);
        }
    }
    cout << "L'algorithme a utilisé " << compteur << " itération(s)" << endl;

    compteur = 0;
    return m;
}

// algo 3 : algo de Wiedemann
Vecteur Wiedemann(Matrice *A, Vecteur b){
    Poly m = pol_min_matrice(A, b);
    int i;
    
    Poly h, h1, h2, x, m0, moins1, tmp;
    m0=creer_poly_nul(b.corps);
    x=creer_poly_nul(b.corps);
    moins1=creer_poly_nul(b.corps);
    
    set(&moins1, 0, b.corps-1);
    set(&x, 1, 1);
    
    set(&m0, 0, get(&m, 0)); // m(0)
    
    h1=diff_poly(&m, &m0);    
    h2 = mult_poly(&m0, &x);
    h=div_poly(&h1, &h2, &tmp);
    h=mult_poly(&h, &moins1);
    
    Vecteur hm=horner_mat(h, A, b);
    
    return hm;    
}

Poly vecteur_to_poly(Vecteur *v){
    Poly p=creer_poly_nul(v->corps);
    int i;
    
    for(i=0; i<v->n; i++)
        set(&p, i, get(v, i));
    
    return p;
}



// ------------------------------------------------------------------//
//                                 MAIN
// ------------------------------------------------------------------//

int main(int argc, char** argv) {

    unsigned int tm[][3]={{0,0,1},{0,1,2},{0,2,3},{1,0,4},{1,2,1},{2,0,1},{2,1,3},{2,2,1}};
    
    Matrice M=creer_matrice(tm, 8, 3, 5);
    
    print_matrice(&M);
    
    unsigned int tb[][2]={{1, 1}, {2, 2}};
    Vecteur b=creer_vecteur(tb, 2, 3, 5);
    
    print_vecteur(&b);
    
    Vecteur res = Wiedemann(&M, b);
    cout << "Résultat : " << endl;
    print_vecteur(&res);
    
    /*Matrice M = creer_matrice_nulle (TAILLE, CORPS);
    while (M.nb_coefficients <= 100/TAUXREMPLISSAGE) {
        int i = rand()%TAILLE;
        int j = rand()%TAILLE;
        if (get(&M, i, j)==0) {
            int v = rand()%CORPS;
            set2 (&M, i, j, v);
        }
    }
    
    Vecteur b = creer_vecteur_nul(TAILLE, CORPS);
    for (int i = 0; i < TAILLE; i++) {
        int v = rand()%CORPS;
        set(&b, i, v);
    }
    
    print_matrice(&M);
    print_vecteur(&b);
    
    Vecteur r=Wiedemann(&M, b);
    
    
    print_vecteur(&r);*/
    
    return 0;
}
