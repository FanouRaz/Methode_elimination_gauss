#include<iostream>
#include<vector>
#include<cmath>
using namespace std;

//a
float const epsilon = pow(10,-6); 

void afficherLigne(vector<float> tab){
    for(float elt : tab){
        cout<<elt<<"\t";
    }
    cout<<endl;
}
//Afficher une matrice n x n
void afficher(vector<vector<float> > matrice){
    for(int i=0;i<matrice.size();i++){
        for(float ligne: matrice[i]){
            cout << ligne << "\t\t";
        }
        cout<< endl << endl;
    }
}

//Fonction pour rechercher la ligne à considerer comme pivot, contenant la valeur maximal sur la colonne col
int positionMaxColonne(vector<vector<float>> tab, int col){
    int pos(-1),max(0);
    for(int i(col);i<tab.size();i++){
        if(tab[i][col] > max){
            max = tab[i][col];
            pos = i;
        }
    }
    return pos;
}

//Implémentation de l'algorithme de l'élimination de Gauss
void gauss(vector<vector<float>> A, vector<float> B){
    //L'algorithme retourne la solution x d'un système d'équation mise sous forme matricielle Ax = b
    vector<float>* Ppivot = nullptr; // Pointeur pointant vers la ligne qui sera le pivot
    vector<float>* Pmax = nullptr;// Pointeur pour récuperer la ligne contenant l'eventuel pivot
    
    /*
        Fonctionnement de l'algorithme:
        Pour i allant de la premiere ligne à l'avant-dernière ligne de A:
            *Trouver l'indice k tel que max(A[k][i]) pour k allant de la ligne i à la dernière ligne de A,
            *Si k != i on interchange la ligne k et la ligne i
            *pour j allant de la ligne i+1 à la dernière ligne:
                on remplace la ligne A[j] par A[j] - A[j][i]*A[i]/A[i][i]

        Remarquons que les opérations élémentaires sur A affectent également la matrice B
    */

    for(int i(0);i<A.size() -1;i++){
        vector<float> temp({});
        float tempB(0);

        Pmax = &A[positionMaxColonne(A,i)];
        
        //Permutation
        if(Pmax != &A[i]){
            //Permutation du tableau A
            temp = A[i];
            A[i] = *Pmax;
            *Pmax = temp;

            //Permutation du tableau B
            tempB = B[positionMaxColonne(A,i)];
            B[positionMaxColonne(A,i)] = B[i];
            B[i] = tempB;

        }

        Ppivot = &A[i];  
        
        for(vector<float>* Piter(&A[i+1]);Piter != &A[A.size()]; Piter++){
            //On stocke (*Piter)[i] dans une variable car lors de l'itération elle sera annuler dès que k = i

            float coef((*Piter)[i]);    

            for(int k(0);k<A.size();k++){
                (*Piter)[k] -=  (*Ppivot)[k]  * coef/ (*Ppivot)[i];

                //Modification de B
                if(k == i){    
                    cout<< B[k] <<" - "<< B[k] << " * " << coef <<"/"<<(*Ppivot)[i]<<" = "<<  B[k] -  B[i]  * coef/ (*Ppivot)[i] <<endl;
                    B[k] -=  B[i]  * coef/ (*Ppivot)[i];        
                }
            }

            
            
        }
    }
    afficher(A);
    cout<<"B:"<<endl;
    afficherLigne(B);
}

int main(){
    
    vector<vector<float>> tab = {{3,8,13},{4,8,12},{2,9,18}};
    vector<vector<float>> tab1 = {{1,2,3},{2,4,5},{1,8,0}};
    vector<float> b({4,5,11});
    afficher(tab1);
    cout<<endl<<endl;
    gauss(tab,b);

    return 0;
}