#include <cstdlib>						// C standard library
#include <cstdio>						// C I/O (for sscanf)
#include <cstring>						// string manipulation	
#include <fstream>					// file I/O
#include <ANN/ANN.h>				// ANN declarations
#include <sstream> 

using namespace std;					// make std:: accessible

int main(){
	// NUMBERS etc...
	cout << "Hello world" << endl;
	const double PI = 3.1416;
	char myGrade='A'; // 1 byte
	bool isHappy=true;  // make them start with "is"!!
	int myAge=39;  // 4 bytes
	float favNum=3.1416; // have digits (up to 8)
	double otherfavNum=1.6180339887; // digits up to 16
	int largestInt=2147483647; // next one is -2147...
	number+=5 // number=number+5
	// && || !!  
	int greetingOption=2;
	switch(greetingOption){
		case 1:
		cout<<"bonjour"<<endl;
		break;
		case 2:
		cout<<"hola"<<endl;
		break;
		default
		cout<<"hello"<<endl;
	}
	int largestnumber=(5<2)?5:2;
	int badNums[5]={4,13,14,24,1}; // cannot change the amount of numbers you can store later.
	cout << "Bad Number 1:"<<badNums[0]<<endl;
	char myName[3][3]={{'L','i','o'}{'T','r','e'}; 

	// while instead of for loops when you don't know ahead of time when loop is going to end
	int randNum=(rand()%100)+1;
	while(randNum!=100){
		cout<<randNum<<", ";
		randNum=(rand()%100)+1;
	}
	// other way of stating while
	do{
		cout<<"Guess btw 1 and 10";
		getline(cin,numberGuessed);
		intNumberGuessed=stoi(numberGuessed);
		cout<<intNumberGuessed<<endl;
	}
	while(intNumberGuessed!=4);

	vector <int> lotteryNumVec(10);
	int lotteryNumVec[5]={4,13,14,24,34};
	lotteryNumVec.insert(lotteryNumVec.begin(),lotteryNumArray,lotteryNumArray+3)

	 int getFactorial(int number){
	 	int sum;
	 	if(number==1)sum=1;
	 	else sum=getFactorial(number-1)*number;
	 	return sum;
return sum
	 }

}
