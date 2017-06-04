#include <cstdlib>
#include <iostream>
#include <tgmath.h>

using namespace std;

struct Poly {
	unsigned int* p; // contient les coefficients de p
        unsigned int corps; // corps du polynôme p
        unsigned int degre; // degré du polynôme p
        unsigned int taille; // nombre de bits de p;
        unsigned int taille_bloc; // nombre de bits pour un coefficient
        unsigned coefficients_par_case; // nombre de coefficients par case
};

/*--------------------------------------------------------------------*/
/*Fonctions sur les polynômes*/
/*--------------------------------------------------------------------*/


// récupère de le coefficient de la puissance 10^i-ème dans le polynôme a
int get(Poly *p, unsigned int i){
    // masque permemttant de récupérer le coefficient :
    unsigned int masque=0;
    int j=0;
    
    /* le masque contient autant de 1 que de bits 
     * réservés à chaque coefficient : */
    
    do{
        masque<<=1;
        masque|=1;
    }while((++j)<p->taille_bloc);
    
    // numéro du bloc où se trouve le coefficient cherché :
    int num_bloc=i/p->coefficients_par_case;
    
    // on récupère le block où se trouve le coefficient :
    unsigned int bloc=p->p[num_bloc];
    
    // on renvoie le coefficient : 
    return (bloc>>((i*p->taille_bloc)%(p->coefficients_par_case*p->taille_bloc)))&masque;
}

int set(Poly *p, unsigned int i, unsigned int a){
    int j=0;
    int b, num_bit;
    
    int num_bloc=i/p->coefficients_par_case;
    int num_coeff=i-num_bloc*p->coefficients_par_case;
    
    num_bit=num_coeff*p->taille_bloc;
    
    while(j<p->taille_bloc){
        b=(a&1);
        if (b)
            p->p[num_bloc] = p->p[num_bloc] | (1u<<(num_bit+j));
        else
            p->p[num_bloc]=p->p[num_bloc] & (~(1u<<(num_bit+j)));
        j++;
        a>>=1;
    }    
}

// ça va pas : 
void print_poly (Poly *p){
    int j=0, a;
    
    /* on affiche le coefficient pour le degré 0 :*/
    // s'il existe et si le coefficient n'est pas 0
    if((j++)<=p->degre && (a=get(p,0))){ 
        cout << a;
    }
    
    /* on affiche le coefficient pour le degré 1 :*/
    if((j++)<=p->degre){
        if(a)
            cout << "+";
        if(a=get(p,1)){
            if (a!=1)
                cout << a;
            cout << "X";
        }
    }
    
    /* on affiche le reste des coefficients :*/
    while(j<=p->degre){
        if(a)
            cout << "+";
        if(a=get(p,j)){
            if(a!=1)
                cout << a;
            cout << "X^" << j;
        }
        
        j++;
    }    
}

Poly creer_poly(unsigned int *t, int degre, int corps){
    int i=0, j=0;
    unsigned int blocs;
    Poly p;
    p.corps=corps;
    p.degre=degre;
    
    /* On arrondit au-dessus pour trouver le nombre nécessaire
     * de bits pour un coefficient :*/
    p.taille_bloc=ceil(log2(corps)); 
    
    // Le nombre de coefficient qu'on peut mettre par case du tableau :
    p.coefficients_par_case = (sizeof(unsigned int))*8/p.taille_bloc;
    
    blocs=degre/p.coefficients_par_case;
    
    // Le nombre de bits :
    p.taille=blocs*sizeof(unsigned int)*8;
    
    
    //On remplit le tableau avec les coefficients du polynôme (bit-packing) :
    
    p.p = (unsigned int*)malloc(p.taille);
   
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

int main(int argc, char** argv) {
    Poly p;
    
    unsigned int t[]={2, 4, 3}, t2[]={3}, t3[]={0, 2, 4, 2, 0, 1}
    , t4[]={2, 3, 4, 1 ,0 , 0 , 1, 0, 4, 7, 5, 8, 6, 3, 6, 7, 12, 8, 3, 4, 4, 1};
    
    p=creer_poly(t4, 21, 23);
    
    print_poly(&p);
    return 0;
}

