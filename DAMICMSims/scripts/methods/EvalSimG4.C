#include "TArrayD.h"
#include "TArrayI.h"
#include "TList.h"
#include "TMath.h"
#include <math.h>

Int_t InsideOutside(vector<Int_t>* Number, TArrayI* in, Double_t  Num){
  if (Number->size() >1){
    if(Num == 10){
      return 2;
    }
    else{
      return in->GetAt(Number->at(0));
    }
  }
  else if(Number->size() == 1){
    return in->GetAt(Number->at(0));
  }
  else {
    return 3;
  }
}

Double_t Process (Int_t number1, Int_t number2){
  if (number1 == -1){
    return 0;
  }
  else if(number1 == 2){
    switch(number2) {
    case 2:
      return 1;
    case 3:
      return 2;
    case 10:
      return 3;
    case 11:
      return 4;
    case 12:
      return 5;
    case 13:
      return 6;
    case 14:
      return 7;
    }
    return 9;
  }
  else if(number1 == 6){
    return 8;
  }
  else{
    return 9;
  }
}

Double_t ClusterProcessOrigin(TArrayI* Num1, TArrayI* Num2, vector<Int_t>* Number){

  if (Number->size() >1){
    Double_t Principal = Process(Num1->GetAt(Number->at(0)), (Int_t)Num2->GetAt(Number->at(0)));

    Bool_t Other = false;
    for (unsigned int i = 1 ; i < Number->size(); i++){
      Double_t New = Process(Num1->GetAt(Number->at(i)), (Int_t)Num2->GetAt(Number->at(i)));
      if( New != Principal && New != 5 && New != 1){
        Other = true;
      }
      if ( New == 8){
        Principal = 8;
      }
    }
    if (Other == true){
      if (Principal != 8){
        return 10;
      }
      else{
        return Principal;
      }
    }
    else{
      return Principal;
    }
  }
  else if (Number->size() == 1){
    return Process(Num1->GetAt(Number->at(0)), (Int_t)Num2->GetAt(Number->at(0)));
  }
  else{
    return 9;
  }
}

Bool_t Match(Double_t Pixel1, Double_t Pixel2){
  Double_t Diff = TMath::Abs(Pixel2 -Pixel1);
  if (Diff>=4){
    return false;
  }
  else{
    return true;
  }
}

vector<Int_t>* EvalSimG4Num(TArrayD* Xas, TArrayD* Yas, TArrayD* Xpart,TArrayD* Ypart){
  Int_t sizeXas = Xas->GetSize();
  Int_t sizeYas = Yas->GetSize();
  Int_t sizeXpart = Xpart->GetSize();
  Int_t sizeYpart = Ypart->GetSize();
  vector<Int_t>* Num = new vector<Int_t>();

  for (Int_t i =0; i < sizeXpart; i++){
    bool found = false;
    int j =0;
    while( j < sizeXas && found == false){
      Bool_t MatchX = Match(Xas->At(j), Xpart->At(i));
      Bool_t MatchY = Match(Yas->At(j), Ypart->At(i));
      if (MatchX && MatchY){
        Num->push_back(i);
        found = true;
      }
      j++;
    }
  }

  return Num;
}


void EvalSimG4Prop(TList* vals, TArrayI* p1, TArrayI* p2, TArrayI* in, vector<Int_t>* Number){
  Double_t Process = ClusterProcessOrigin(p1, p2, Number);
  AddVariable(vals,"Process",Process);
  Double_t InsideOut = (Int_t)InsideOutside(Number, in, Process);
  AddVariable(vals, "Inside", InsideOut);
}
