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
    
    // degré de res
    unsigned int d=max(a->degre, b->degre);
    
    // tableau de coefficients nuls
    unsigned int *t=(unsigned int*)malloc(max(a->taille, b->taille)); 
    unsigned int blocs = max(a->blocs, b->blocs);
    int i;
    for(i=0; i<=blocs; i++)
        t[i]=0;
    
    res=creer_poly(t, d, a->corps);
    
    res.corps=a->corps;
    res.degre=d;
    
    for(i=res.degre; i>=0; i--)
        set(&res, i, ((get(a, i)+res.corps-get(b, i))%res.corps)); 
    
    res.coefficients_par_case=a->coefficients_par_case;
    res.taille_bloc=a->taille_bloc;
    
    // on met à jour res :
    unsigned int k=0;
    while(res.degre>=0 && get(&res, res.degre)==0){
        res.degre--;
        k++;
    }
    res.blocs=blocs-k/res.taille_bloc;
    res.taille=8*sizeof(unsigned int);
            
    res.p=(unsigned int*)realloc(res.p, res.taille);
    
    return res;
}

// Renvoie l'addition de deux polynômes à valeurs dans le même corps
Poly add_poly(Poly *a, Poly *b){
    Poly res;
    
    // degré de res
    unsigned int d=max(a->degre, b->degre);
    
    // tableau de coefficients nuls
    unsigned int *t=(unsigned int*)malloc(max(a->taille, b->taille));
    unsigned int blocs = max(a->blocs, b->blocs);
    int i;
    for(i=0; i<=blocs; i++)
        t[i]=0;
    
    res=creer_poly(t, d, a->corps);
    
    res.corps=a->corps;
    res.degre=d;
    
    for(i=res.degre; i>=0; i--)
        set(&res, i, ((get(a, i)+get(b, i))%res.corps));
    
    res.coefficients_par_case=a->coefficients_par_case;
    res.taille_bloc=a->taille_bloc;
    
    // on met à jour res :
    unsigned int k=0;
    while(res.degre>=0 && get(&res, res.degre)==0){
        res.degre--;
        k++;
    }
    res.blocs=blocs-k/res.taille_bloc;
    res.taille=8*sizeof(unsigned int);
    
    res.p=(unsigned int*)realloc(res.p, res.taille);
    
    return res;
}



//Multiplication de deux polynômes à coefficients dans le même corps
Poly mult_poly(Poly *a, Poly *b){
    Poly res;
    
    // degré de res
    unsigned int d=(a->degre)*(b->degre);
    
    // tableau de coefficients nuls
    unsigned int *t=(unsigned int*)malloc((a->taille)*(b->taille));
    unsigned int blocs = (a->blocs)*(b->blocs);
    int i;
    for(i=0; i<=blocs; i++){
        t[i]=0;
    }
    
    res=creer_poly(t, d, a->corps);
    
    res.corps=a->corps;
    res.degre=d;
    
    //set(&res, d, (get(a, a->degre)+1)*(get(b, (b->degre))+1)); // coeff de plus haut degré
    for(int i=0; i<=a->degre; i++){
        
        for (int j=0; j<=b->degre; j++) {
            int r = get(&res, i+j); // valeur du coeff à la puissance 1O^{i+j}
            if (r==0) {
                set(&res, i+j, (get(a, i)*get(b, j))%res.corps);
            }
            else {
                set(&res, i+j, (r+ get(a, i)*get(b, j))%res.corps);
            }
        }
    }
    
    res.coefficients_par_case=a->coefficients_par_case;
    res.taille_bloc=a->taille_bloc;
    
    // on met à jour res :
    unsigned int k=0;
    while(res.degre>=0 && get(&res, res.degre)==0){
        res.degre--;
        k++;
    }
    res.blocs=blocs-k/res.taille_bloc;
    res.taille=8*sizeof(unsigned int);
    
    res.p=(unsigned int*)realloc(res.p, res.taille);
    
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
    
    return 0;
}

