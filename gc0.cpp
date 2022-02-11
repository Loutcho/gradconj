/*
METHODES PERFORMANTES EN CALCUL SCIENTIFIQUE
Méthode du gradient conjugué
appliquée au problème :
	-laplacien(u)=f sur U=]0,1[x]0,1[
   u=0 sur frontière(U)
Luc ROUSSEAU
30.10.1997
*/

#include <math.h>  // appel à sinus
#include <stdio.h> // appel à fprintf

const float pi=3.14159265359;

const int n=100;				  // ne pas dépasser 125, sinon mémoire insuffisante
const int N=10000;		     // N=n*n
const float h=1/101.0;       // h=1/(n+1)
const float epsilon=1/100.0; // epsilon ~= h

class vecteur
{
private:
	float tab[N];
public:
	// Initialisation du second membre b de l'équation Ax=b à résoudre,
	// selon la syntaxe b.secondmembre(ff);
   // Comme effet de bord, ff est le vecteur représentant la fonction f
   void secondmembre(vecteur&);

   // Initialisation d'un vecteur à zéro
   void initialisation();

	// Somme vectorielle
	vecteur operator+(vecteur);

	// Produit externe d'un vecteur par un scalaire
   friend vecteur operator*(float,vecteur);

   // Produit scalaire entre deux vecteurs
   float operator*(vecteur);

   // Multiplication par A d'un vecteur
   friend vecteur A(vecteur);

   // Affichage d'un vecteur (dans un fichier)
   void affichage(FILE *);

   // Rapport de deux vecteurs
   // Utile pour tester si ff/u = 2*pi*pi (cf. choix de f)
   vecteur operator/(vecteur);
};

float f(float X,float Y)
{
	return 1000.0*sin(pi*X)*sin(pi*Y);
   // Cette fonction f est fonction propre de l'opérateur laplacien.
   // Comme u=f/(2*pi*pi) vérifie -laplacien(u)=f et u(bords)=0,
   // un test de validité de l'algorithme est : a-t-on f/u=2*pi*pi=19.739 ?
}

void vecteur::secondmembre(vecteur &ff)
{
	int i,voisin;
   float xx,yy,q,hhq;
   for (i=0;i<N;i++)
   {
   	xx=h+(i%n)*h;
      yy=h+(i/n)*h;
      q=f(xx,yy);
      hhq=h*h*q;
      voisin=0;
      if (i%n !=  0  )  { voisin++; };
      if (i-n >=  0  )  { voisin++; };
      if (i%n != n-1 )  { voisin++; };
      if (i+n  <  N  )  { voisin++; };
      switch (voisin)
      {
      case 4:tab[i]=hhq; break;
      case 3:tab[i]=2.0*hhq/3.0; break;
      case 2:tab[i]=hhq/3.0; break;
      }
      ff.tab[i]=q;
   }
}

void vecteur::initialisation()
{
	int i;
   for (i=0;i<N;i++)
   {
   	tab[i]=0.0;
   }
}

vecteur vecteur::operator+(vecteur v)
{
	vecteur vv;
	int i;
   for (i=0;i<N;i++)
   {
   	vv.tab[i]=tab[i]+v.tab[i];
   }
   return vv;
}

vecteur operator*(float a,vecteur v)
{
	vecteur vv;
   int i;
   for (i=0;i<N;i++)
   {
   	vv.tab[i]=a*v.tab[i];
   }
   return vv;
}

float vecteur::operator*(vecteur v)
{
	float P=0.0;
   int i;
   for (i=0;i<N;i++)
   {
   	P+=tab[i]*v.tab[i];
   }
	return P;
}

vecteur A(vecteur v)
{
	vecteur vv;
   int i;
   for (i=0;i<N;i++)
   {
   	vv.tab[i]=4*v.tab[i];
      if (i%n !=  0  )  { vv.tab[i] -= v.tab[i-1]; };
      if (i-n >=  0  )  { vv.tab[i] -= v.tab[i-n]; };
      if (i%n != n-1 )  { vv.tab[i] -= v.tab[i+1]; };
      if (i+n  <  N  )  { vv.tab[i] -= v.tab[i+n]; };
   }
   return vv;
}

void vecteur::affichage(FILE *stream)
{
	int ligne,colonne;
   fprintf(stream,"{\n");
   for (ligne=n-1;ligne>=0;ligne--)
   {
   	for (colonne=0;colonne<n;colonne++)
      {
   		fprintf(stream,"%6.3f ",tab[colonne+n*ligne]);
      }
      fprintf(stream,"\n");
   }
   fprintf(stream,"}\n");
}

vecteur vecteur::operator/(vecteur v)
{
	vecteur vv;
   int i;
   for (i=0;i<N;i++)
   {
   	vv.tab[i]=tab[i]/v.tab[i];
   }
   return vv;
}

void main(void)
{
	FILE *stream;
   stream = fopen("c:\\windows\\bureau\\gc0.txt","w+");

	vecteur b,x,r,p,z,ff;
   float alpha,beta,delta0,delta1;
   int iteration=0;
   b.secondmembre(ff);
 	fprintf(stream,"fonction f de départ, discrétisée : \n");
   	ff.affichage(stream);
   x.initialisation();

   r=A(x)+(-1)*b;
 	p=(-1)*r;
   delta0=r*r;

   fprintf(stream,"pas | delta0/epsilon\n");
   fprintf(stream,"----+---------------\n");

   while (delta0>epsilon)
   {
   	iteration++;
   	z=A(p);
      alpha=delta0/(z*p);
      x=x+alpha*p;
      r=r+alpha*z;
      delta1=r*r;
      beta=delta1/delta0;
      p=(-1)*r+beta*p;
      delta0=delta1;
      fprintf(stream,"%03d | %f\n",iteration,delta0/epsilon);
   }
   fprintf(stream,"%d iterations\n",iteration);
   fprintf(stream,"solution x : \n"); x.affichage(stream);
	fprintf(stream,"test de validité : rapport ff/x : \n");
   	(ff/x).affichage(stream);
   fclose(stream);
}
