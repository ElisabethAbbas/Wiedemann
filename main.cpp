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
        unsigned int blocs; // nombre de blocs au total
};

/*--------------------------------------------------------------------*/
/*Fonctions sur les polynômes*/
/*--------------------------------------------------------------------*/


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
    }while((++j)<p->taille_bloc);
    
    // numéro du bloc où se trouve le coefficient cherché :
    int num_bloc=i/p->coefficients_par_case;
    
    // on récupère le block où se trouve le coefficient :
    unsigned int bloc=p->p[num_bloc];
    
    // on renvoie le coefficient : 
    return (bloc>>((i*p->taille_bloc)%(p->coefficients_par_case*p->taille_bloc)))&masque;
}

// modifie le coefficient de X^i dans le polynôme *p par a
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

/* crée et renvoie le polynôme p via le tableau de coefficients "t", 
 * de degré "degre" et sur le corps "corps". */
Poly creer_poly(unsigned int *t, int degre, int corps){
    int i=0, j=0;
    Poly p;
    p.corps=corps;
    p.degre=degre;
    
    /* On arrondit au-dessus pour trouver le nombre nécessaire
     * de bits pour un coefficient :*/
    p.taille_bloc=ceil(log2(corps)); 
    
    // Le nombre de coefficient qu'on peut mettre par case du tableau :
    p.coefficients_par_case = (sizeof(unsigned int))*8/p.taille_bloc;
    
    p.blocs=degre/p.coefficients_par_case;
    
    // Le nombre de bits :
    p.taille=p.blocs*sizeof(unsigned int)*8;
    
    
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

// renvoie la différence du polynôme a par le polynôme b, 
// il faut que les deux polynômes aient des coefficients dans le même corps
Poly diff_poly(Poly *a, Poly *b){
    Poly res;
    unsigned int t[3]={0, 0, 0};
    int i;
    
    res=creer_poly(t, max(a->degre, b->degre), a->corps);
    
    for(i=0; i<=res.degre; i++){
        set(&res, i, (get(a, i)+res.corps-get(b, i))%res.corps); 
    }
    
    return res;
}

// Renvoie l'addition de deux polynômes à valeurs dans le même corps
Poly add_poly(Poly *a, Poly *b){
    Poly res;
    int d = max(a->degre, b->degre);
    unsigned int t[d+1];
    for (int i = 0; i<=d+1; i++) {
        t[i] = 0;
    }
    res = creer_poly(t, d, a->corps);
    
    for (int i = 0; i <=(res.degre); i++) {
        set(&res, i, (get(a,i)+get(b,i))%(res.corps));
    }
    return res;
}

// Trouver le polynôme réciproque revP(x) = x^d P(1/x)
Poly rev(Poly p){
    unsigned int t[p.degre];
    for (int i = 0; i<=(p.degre); i++) {
        int a = get(&p,i);
        t[(p.degre)-i]=a;
    }
    Poly rev_p = creer_poly(t, p.degre, p.corps);
    return rev_p;
}


// ------------------------------------------------------------------//
//                                 MAIN
// ------------------------------------------------------------------//

int main(int argc, char** argv) {
    Poly p, p1, p2, p3, p4, p5, p6;
    
    unsigned int t[]={2, 4, 3};
    unsigned int t2[]={3};
    unsigned int t4[]={2, 3, 4, 1 ,0 , 0 , 1, 0, 4, 7, 5, 8, 6, 3, 6, 7, 12, 8, 3, 4, 4, 1};
    unsigned int t5[]={1, 4, 6};

    
    p=creer_poly(t4, 21, 23);
    p1=creer_poly(t, 2, 7);
    p2=creer_poly(t5, 2, 7);
    p3=diff_poly(&p1, &p2);
    

    cout << "P(x) = ";
    print_poly(&p);
    Poly rev_p = rev(p);
    cout << endl << "Polynôme réciproque de P : "<< endl;
    print_poly(&rev_p);    
    
    cout << endl;
    
    cout << "Le troisième polynôme est le résultat de la différence des deux premiers (corps Z/7Z) :" << endl;
    print_poly(&p1);
    print_poly(&p2);
    print_poly(&p3);
    
    
    // TEST ADDITION:
    unsigned int t3[]={1,2,4,2,0,1};
    unsigned int t6[]={1,3,4,1,0,2,1,1,3};
    p4=creer_poly(t3, 5, 5);
    p5=creer_poly(t6, 8, 5);
    p6=add_poly(&p4, &p5);
    cout << endl << "ADDITION DANS Z/5Z : " << endl;
    print_poly(&p4);
    cout << "+ " << endl;
    print_poly(&p5);
    cout << "= " << endl;
    print_poly(&p6);
    cout << endl;
    
    return 0;
}

