/*
    RAZAKAHASINA Fanomezana Sarobidy
    L3 MISA 
    Méthode d'élimination de Gauss
*/

#include<iostream>
#include<vector>
#include<cmath>
#include<fstream>
using namespace std;

int const EPSILON = pow(10, -6);


vector<float> remplirLigne(string ligne, int n);
void lireData(int &n,vector<vector<float>> &A, vector <float> &b);
void gauss(vector<vector<float>> A, vector<float> B);
void afficherVect(vector<float> tab);
int positionMaxColonne(vector<vector<float>> tab, int col);
void afficher(vector<vector<float> > matrice);

int main(){
    int n = 0;
    vector<vector<float>> A({});
    vector<float> B({});
    lireData(n,A,B);
    cout<<"Ce programme permet de résoudre un système d'equation par la méthode de l'élimination de Gauss en le mettant"
    <<"sous forme d'équation matricielle Ax = B. L'utilisateur inserera les matrices A et B dans le fichier data.txt et le programme"
    <<"affichera les solutions"<<endl;
    cout<<"La matrice A proposée par l'utilisateur est: "<<endl;
    afficher(A);
    cout<<"La matrice B proposée par l'utilisateur est:"<<endl;
    afficherVect(B);
     
    cout<<endl<<endl;
    gauss(A,B);

    return 0;
}

//Récuperation des matrices dans le fichier data.txt
void lireData(int &n,vector<vector<float>> &A, vector <float> &b){
    fstream donnee;    
    string line;
    int taille(0);
    int tailleA=0;
    int tailleB=0;
    //ouverture du fichier en mode lecture
    donnee.open("data.txt", ios::in);                                
    if(!donnee){
        cout<<"Erreur lors du chargement des données. Veuillez placer le fichier data.txt dans le répertoire courant."<<endl;
        exit(-1);                                                                                       
    }
    
    while(getline(donnee,line)){
        if(taille==0)  n=stoi("0"+line);
        else if(tailleA < n){
            A.push_back(remplirLigne(line , n));
            tailleA++;
        }
        else{
            b.push_back(stof(line));
            tailleB++;
        }
        taille++;
    }

    donnee.close();
    if(tailleA!= n){
        cout<<"La dimension de la matrice A dans le fichier 'data.txt' n'est pas valide"<<endl;
        exit(-1);
    }
    if(tailleA != tailleB){
        cout<<"La dimension du vecteur B dans le fichier 'data.txt' n'est pas valide"<<endl;
        exit(-1);
        }
    }

//Construire les matrices A et b à partir des données recupérées dans le fichier data.txt 
vector<float> remplirLigne(string ligne, int n){
    vector <float> row;
    int i=0;
    int colonne=0;
    ligne+=" 0";
    int len =ligne.length();
    string valeurString="";
    while(i<len){
        if((ligne.at(i)>='0'&& ligne.at(i)<='9')||(ligne.at(i)=='.')||(ligne.at(i)=='+')||(ligne.at(i)=='-')){
            valeurString+=ligne.at(i);
        }
        else{
            if(valeurString.length()){
                row.push_back(stof(valeurString));
                colonne++;
            }
            valeurString="";
        }
            i++;
    }
    if(colonne != n){
        cout<<"La dimension de la matrice dans le fichier 'data.txt' n'est pas valide"<<endl;
        exit(-1);
    }
    return row;
}

//Fonction pour afficher un vecteur colonne
void afficherVect(vector<float> tab){
    for(float elt : tab){
        cout<<elt<<endl;
    }
    cout<<endl;
}

//Afficher une matrice n x n
void afficher(vector<vector<float> > matrice){
    for(unsigned int i=0;i<matrice.size();i++){
        for(float ligne: matrice[i]){
            cout << ligne << "\t\t";
        }
        cout<< endl << endl;
    }
}


//Solution x de l'equation Ax = B, A étant une matrice inversible et B un vecteur colonne 
void solution (vector<vector<float>> A,vector<float> B){
    float x[A.size()] = {0};
    float somme;
    for(int i(A.size() -1);i>=0;i--){
        somme = 0;
        for(int j(i+1);j<A.size();j++){
            somme += (x[j] * A[i][j]);
        }
        x[i]=((B[i] -somme)/ A[i][i]);
    }
    
    cout << "L'équation admet pour solution : (";
    for(int j(0);j<A.size();j++){
        cout<< x[j];
        if(j < A.size() - 1){
            cout << " , ";
        } 
    }
    cout<<")"<<endl;
}

//Fonction pour rechercher la ligne à considerer comme pivot, contenant la valeur maximal sur la colonne col
int positionMaxColonne(vector<vector<float>> tab, int col){
    int pos(-1),max(0);
    for(int i(col);i<tab.size();i++){
        if(fabs(tab[i][col]) > max){
            max = fabs(tab[i][col]);
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
        Fonctionnem4.17233e-07ent de l'algorithme:
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
            //Permutation du tableau B
            tempB = B[positionMaxColonne(A,i)];
            B[positionMaxColonne(A,i)] = B[i];
            B[i] = tempB;

            //Permutation du tableau A
            temp = A[i];
            A[i] = *Pmax;
            *Pmax = temp;
            
        
        }

        Ppivot = &A[i];  
        int c(i+1);
        for(vector<float>* Piter(&A[c]);Piter != &A[A.size()]; Piter++){
            //On stocke (*Piter)[i] dans une variable car lors de l'itération elle sera annuler dès que k = i
            float coef((*Piter)[i]/(*Ppivot)[i]);    
        
            //Modification de B
            B[c] -=  B[i]  * (*Piter)[i]/ (*Ppivot)[i];

            //Modification de A
            for(int k(0);k<A.size();k++){
                (*Piter)[k] -=  (*Ppivot)[k]  * coef;
            }
            
            c++;
        }

    }
    cout<<"La matrice A une fois l'algorithme achevée:"<<endl;
    afficher(A);
    cout<<"La matrice B une fois l'algorithme achevée:"<<endl;
    afficherVect(B);
    solution(A,B);
}
