#include "header.h"


// check the validity of the parameters' range:
double CheckCor(int *cormod, double *par)
{
  double col=0.0, R_power=0.0, R_power1=0.0, R_power2=0.0, R_power_s=0.0, R_power_t=0.0, var11=0.0, var22=0.0;
  double rho=0.0, sep=0, scale=0.0, smooth=0.0, scale_s=0.0, scale_t=0,smooth_s=0.0,smooth_t=0.0;
  double scale11=0.0, scale22=0, scale12=0.0, smoo11=0.0, smoo12=0.0, smoo22=0,
    R_power11=0.0, R_power12=0.0, R_power22=0,nug11=0.0,nug22=0.0;


  switch(*cormod) // Correlation functions are in alphabetical order
    {
    case 1:// Cauchy correlation function
      R_power2=par[0];
      scale=par[1];
      if(scale<=0 || R_power2<=0) rho=-2;
     break;
    case 2:
    case 3:
    case 4:// Exponential correlation function
    //case 6:// Gaussian correlation function
    case 10:// skarofsky
    case 16://wave correlation function
      scale=par[0];
      if(scale<=0) rho=-2;
      break;
    case 8: // Generalised Cuachy correlation function
    case 5: // Dagum
      R_power1=par[0];
      R_power2=par[1];
      scale=par[2];
      if(scale<=0 || R_power1<=0 || R_power1>2 || R_power2<=0) rho=-2;
      break;
    case 12:// Stable correlation function
    case 17: //multiquadric sphere
      R_power=par[0];
      scale=par[1];
      if(scale<=0 || R_power<0 || R_power>2) rho=-2;
       break;
    case 11:// wen0 correlation function
        R_power=par[0];
        scale=par[1];
        if(scale<=0 || R_power<1.5) rho=-2;
       break;
    case 13:// wen1 correlation function
        R_power=par[0];
        scale=par[1];
        if(scale<=0 || R_power<2.5) rho=-2;
      break;
    case 15:// wen1 correlation function
        R_power=par[0];
        scale=par[1];
        if(scale<=0 || R_power<3.5) rho=-2;
        break;
    case 19: // Generalised wendland
               R_power1=par[0];
        scale=par[1];
        smooth=par[2];
       //if(scale<=0 ||  R_power1>(1.5+smooth) ||smooth<0) rho=-2;
      
            if(scale<=0 ||smooth<-0.5||R_power1>(1.5+smooth) ) rho=-2;
      break;
          case 26: // Generalised wendland hole
          case 29:  
               R_power1=par[0];
        scale=par[1];
        smooth=par[2];
       //if(scale<=0 ||  R_power1>(1.5+smooth) ||smooth<0) rho=-2;
      
            if(scale<=0 ||smooth<-0.5||R_power1>(3.5+smooth) ) rho=-2;
      break;
           case 27: // Matern hole
               R_power1=par[0];
        scale=par[1];
        smooth=par[2];
       //if(scale<=0 ||  R_power1>(1.5+smooth) ||smooth<0) rho=-2;
            if(scale<=0 ||smooth<0.0 ) rho=-2;
      break;
         case 28: // Matern hole
        scale=par[0];
        smooth=par[1];
       //if(scale<=0 ||  R_power1>(1.5+smooth) ||smooth<0) rho=-2;
            if(scale<=0 ||smooth<0.0 ) rho=-2;
      break;
    case 6:
    case 31:    
        R_power1=1/par[0];
        scale=par[1];
        smooth=par[2];
       //if(scale<=0 ||  R_power1>(1.5+smooth) ||smooth<0) rho=-2;
            if(scale<=0 ||smooth<-0.5||R_power1>(1.5+smooth) ) rho=-2;
      break;
      case 24: //kummer
      case 25:  
        R_power1=par[0];
        scale=par[1];
        smooth=par[2];
       //if(scale<=0 ||  R_power1<0 ||smooth<0) rho=-2;
            if(scale<=0 ||smooth<0) rho=-2;
      break;
    case 21: // hyperg   
        R_power=par[0];
         R_power1=par[1];
        scale=par[2];
        smooth=par[3];
       //if(scale<=0 ||smooth<1||((2*(R_power1-smooth)*(R_power-smooth))<smooth)||((2*(R_power1+R_power)<(6*smooth+1)))) rho=-2;
              if(scale<=0 ||smooth<1) rho=-2;
      break;
       case 22:   
       case 23:
       case 30:    
        R_power=par[0];
         R_power1=par[1];
        scale=par[2];
        smooth=par[3];
       //if(scale<=0 ||smooth<1||((2*(R_power1-smooth)*(R_power-smooth))<smooth)||((2*(R_power1+R_power)<(6*smooth+1)))) rho=-2;
              if(scale<=0 ||smooth<= -0.5) rho=-2;
      break;
       case 7:
        R_power1=par[0];
        scale=par[1];
        smooth=par[2];
       if(scale<=0 ||  R_power1<(1.5+smooth)) rho=-2;
            if(scale<=0 ||smooth<0) rho=-2;
      break;
    case 18://sinR_power valid on sphere
            R_power=par[0];
            if( R_power<0 ) rho=-2;
            break;
    case 14://  Whittle-Matern correlation function
    case 20:
      scale=par[0];
      smooth=par[1];
      if(scale<=0 || smooth<=0) rho=-2;
      break;
      // START non-separable correlation functions:
    case 42: //Gneiting correlation model as in (14) Gneitint (2002) with tau=1
    case 46://Porcu model as in (4) of Porcu Bevilacqua et al (2010), with beta_1=beta_2=1
    case 50:
    case 60:
    case 52:
    case 54:
      R_power_s=par[0];
      R_power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      sep=par[4];
      if(scale_s<=0 || scale_t<=0 || R_power_s<0 || R_power_s>2 || R_power_t<0 || R_power_t>2  || sep<0 || sep>1) rho=-2;  //||
      break;
             case 61:
        R_power_s=par[0];
        scale_s=par[2];
        scale_t=par[3];
        R_power=par[1];
        sep=par[4];
        smooth=par[5];
        if(scale_s<=0 || scale_t<=0 || R_power_s<0 || R_power_s>2 ||  sep<0 || sep>1 || smooth<=0||R_power<1) rho=-2;  //||
        break;
        case 62:
         R_power_t=par[0];
        scale_s=par[2];
        scale_t=par[3];
        R_power=par[1];
        sep=par[4];
        smooth=par[5];
        if(scale_t<=0 || scale_s<=0 || R_power_t<0 || R_power_t>2 ||  sep<0 || sep>1 || smooth<=0||R_power<1) rho=-2;  //||
        break;
  case 44:// Iaco-Cesare model as in (14) of Gneitint (2006): note that R_power parameters are in [0,1]
      R_power2=par[0];
      R_power_s=par[1];
      R_power_t=par[2];
      scale_s=par[3];
      scale_t=par[4];
      if(R_power2<=0||scale_s<=0 || scale_t<=0 || R_power_s<0 || R_power_s>2 || R_power_t<0 || R_power_t>2) rho=-2;
      break;
    case 48:// Stein model as in (16) of JASA (2005) with epsilon=0*/
      R_power_t=par[0];
      scale_s=par[1];
      scale_t=par[2];
      smooth=par[3];
      if(scale_s<=0 || scale_t<=0 || R_power_t<0 || R_power_t>2 || smooth<=0) rho=-2;
      break;
    case 58:  //only sphere
    case 56:    //only sphere1
        R_power_s=par[0];
        R_power_t=par[1];
        scale_s=par[2];
        scale_t=par[3];
        if(scale_s<=0  || scale_t<=0 || R_power_s<0 || R_power_s>2 || R_power_t<0 || R_power_s>2) rho=-2;
        break;
    case 63:  // non separable temporal wendloand
    case 65:
    case 67:
        R_power_t=par[0];
        R_power_s=par[1];
        R_power=par[2];
        scale_s=par[3];
        scale_t=par[4];
        sep=par[5];
  if(*cormod==63&&(scale_s<=0||scale_t<=0 || R_power<2.5  || R_power_t>2|| R_power_t<0||sep<0 ||sep >1|| R_power_s<3.5)) rho=-2;
  if(*cormod==65&&(scale_s<=0||scale_t<=0 || R_power<4.5  || R_power_t>2|| R_power_t<0||sep<0 ||sep >1|| R_power_s<4.5)) rho=-2;
  if(*cormod==67&&(scale_s<=0||scale_t<=0 || R_power<6.5  || R_power_t>2|| R_power_t<0||sep<0 ||sep >1||R_power_s<5.5)) rho=-2;
    break;
        case 87:
        R_power_t=par[0];
        R_power_s=par[1];
        R_power=par[2];
        scale_s=par[3];
        scale_t=par[4];
        sep=par[5];
        smooth=par[6];
  if(scale_s<=0||scale_t<=0 || R_power<(2.5+2*smooth)  || R_power_t>2|| R_power_t<0||sep<0 ||sep >1|| R_power_s<(3.5+smooth)||smooth<0) rho=-2;
      break;
   case 64:  // non separable temporal wendloand
    case 66:
    case 68:
       R_power_s=par[0];
        R_power_t=par[2];
        R_power=par[1];
        scale_s=par[3];
        scale_t=par[4];
        sep=par[5];
   if(*cormod==64&&(scale_s<=0||scale_t<=0 || R_power<2.5  || R_power_s>2|| R_power_s<0||sep<0 ||sep >1)) rho=-2;
  if(*cormod==66&&(scale_s<=0||scale_t<=0 || R_power<3.5  || R_power_s>2|| R_power_s<0||sep<0 ||sep >1)) rho=-2;
  if(*cormod==68&&(scale_s<=0||scale_t<=0 || R_power<4.5  || R_power_s>2|| R_power_s<0||sep<0 ||sep >1)) rho=-2;
    break;
        case 88:
       R_power_s=par[0];
        R_power_t=par[2];
        R_power=par[1];
        scale_s=par[3];
        scale_t=par[4];
        sep=par[5];
        smooth=par[6];
   if(scale_s<=0||scale_t<=0 || R_power<(2.5+2*smooth) || R_power_s>2|| R_power_s<0||sep<0 ||sep >1|| R_power_s<(3.5+smooth)||smooth<0)  rho=-2;
    break;
         case 89:
        scale_s=par[0];
        scale_t=par[1];
        smooth_s=par[2];
        smooth_t=par[3];
        sep=par[4];
   if(scale_s<=0||scale_t<=0 || smooth_s<=0||smooth_t<=0 ||  sep<0 )  rho=-2;
    break;

         case 85:
        scale_s=par[0];
        scale_t=par[1];
        R_power_s=par[2];
      R_power_t=par[3];
        smooth_s=par[4];
        smooth_t=par[5];
        sep=par[6];
   if(scale_s<=0||scale_t<=0 || smooth_s<=0||smooth_t<=0 ||  sep<0 || R_power_s<1.5|| R_power_t<1.5)  rho=-2;
    break;

    case 69:  // separable wendlands
    case 70:
    case 71:
    case 72:
    case 73:
    case 74:
    case 75:
    case 76:
    case 77:
          R_power_s=par[0];
          R_power_t=par[1];
          scale_s=par[2];
          scale_t=par[3];
        if(*cormod==69 &&(scale_s<=0  || scale_t<=0  || R_power_s<2.5||R_power_t<2.5)) rho=-2;
        if(*cormod==70 &&(scale_s<=0  || scale_t<=0  || R_power_s<2.5||R_power_t<3.5)) rho=-2;
        if(*cormod==71 &&(scale_s<=0  || scale_t<=0  || R_power_s<2.5||R_power_t<4.5)) rho=-2;
        if(*cormod==72 &&(scale_s<=0  || scale_t<=0  || R_power_s<3.5||R_power_t<2.5)) rho=-2;
        if(*cormod==72 &&(scale_s<=0  || scale_t<=0  || R_power_s<3.5||R_power_t<3.5)) rho=-2;
        if(*cormod==74 &&(scale_s<=0  || scale_t<=0  || R_power_s<3.5||R_power_t<4.5)) rho=-2;
        if(*cormod==75 &&(scale_s<=0  || scale_t<=0  || R_power_s<4.5||R_power_t<2.5)) rho=-2;
        if(*cormod==76 &&(scale_s<=0  || scale_t<=0  || R_power_s<4.5||R_power_t<3.5)) rho=-2;
        if(*cormod==77 &&(scale_s<=0  || scale_t<=0  || R_power_s<4.5||R_power_t<4.5)) rho=-2;
        break;
      // END non-separable correlation functions
      // START separable correlation functions:
    case 82:// Exp-Cauchy:
      R_power2=par[0];
      scale_s=par[1];
      scale_t=par[2];
      if(scale_s<=0 || scale_t<=0 || R_power2<=0) rho=-2;
      break;
 
    case 90:// Matern-Cauchy:
      R_power2=par[0];
      scale_s=par[1];
      scale_t=par[2];
      smooth=par[3];
      if(scale_s<=0 || scale_t<=0 || R_power2<=0 || smooth<=0) rho=-2;
      break;
    case 92:// Matern-exp:
      scale_s=par[0];
      scale_t=par[1];
      smooth=par[2];
      if(scale_s<=0 || scale_t<=0 || smooth<=0) rho=-2;
      break;
    case 94:// Stable-stable:
      R_power_s=par[0];
      R_power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      if(scale_s<=0 || scale_t<=0 || R_power_s<0 || R_power_s>2 || R_power_t<0 || R_power_s>2) rho=-2;
      break;
      // END separable correlation functions:
    case 111:  //wend sep  k=0
    case 113:   //wend sep  k=1
    case 115:   //wend sep  k=2
       var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        R_power=par[5];
        scale=par[6];
     if(*cormod==111&&(nug11<0 || nug22<0 || var11<=0 || var22<=0 || scale<=0 ||
        R_power<2.5||fabs(col)>=1)) { rho=-2;break; }
      if(*cormod==113&&(nug11<0 || nug22<0 || var11<=0 || var22<=0 || scale<=0 ||
        R_power<3.5||fabs(col)>=1)) { rho=-2;break; }
     if(*cormod==115&&(nug11<0 || nug22<0 || var11<=0 || var22<=0 || scale<=0 ||
        R_power<4.5||fabs(col)>=1)) { rho=-2;break; }
     break;
    case 129:  //wend  contr sep  k=0
    case 131:   //wend contr sep  k=1
    case 120:   //wend contr sep  k=2
       var11=par[0];
     var22=par[1];
     nug11=par[2];
     nug22=par[3];
     col=par[4];
     R_power11=par[5];
     R_power22=par[6];
     scale11=par[7];
     scale22=par[8];
      if(*cormod==129&&(nug11<0 || nug22<0 || var11<=0 || var22<=0 || scale11<=0 ||
         scale22<=0 ||  R_power11<2.5||  R_power22<2.5||fabs(col)>= R_pow( R_pow(0.5*(scale11+scale22),2)/(scale11*scale22) ,2.5) )) { rho=-2;break; }
        if(*cormod==131&&(nug11<0 || nug22<0 || var11<=0 || var22<=0 || scale11<=0 ||
        scale22<=0 || R_power11<2.5|| R_power22<2.5||fabs(col)>=1)) { rho=-2;break; }
        if(*cormod==120&&(nug11<0 || nug22<0 || var11<=0 || var22<=0 || scale11<=0 ||
        scale22<=0 ||   R_power11<2.5|| R_power22<2.5||fabs(col)>=1)) { rho=-2;break; }
    break;

    case 112:  //wend full  k=0
    case 114:   //wend full  k=1
    case 116:   //wend full  k=2

    var11=par[0];
    var22=par[1];
    nug11=par[2];
    nug22=par[3];
    col=par[4];
    R_power11=par[5];
    R_power12=par[6];
    R_power22=par[7];
    scale11=par[8];
    scale12=par[9];
    scale22=par[10];
        if(*cormod==112&&(nug11<0 || nug22<0 || var11<=0 || var22<=0 || scale11<=0 ||
        scale12<=0 || scale22<=0 ||
        R_power11<2.5|| R_power12<2.5|| R_power22<2.5||fabs(col)>=1)) { rho=-2;break; }
        if(*cormod==114&&(nug11<0 || nug22<0 || var11<=0 || var22<=0 || scale11<=0 ||
        scale12<=0 || scale22<=0 ||
        R_power11<2.5|| R_power12<2.5|| R_power22<2.5||fabs(col)>=1)) { rho=-2;break; }
        if(*cormod==116&&(nug11<0 || nug22<0 || var11<=0 || var22<=0 || scale11<=0 ||
        scale12<=0 || scale22<=0 ||
        R_power11<2.5|| R_power12<2.5|| R_power22<2.5||fabs(col)>=1)) { rho=-2;break; }
    break;
   case 122:  //bivariate sep matern
     var11=par[0];
     var22=par[1];
     nug11=par[2];
     nug22=par[3];
     col=par[4];
     scale=par[5];
     smooth=par[6];
     if(nug11<0 || nug22<0 || var11<=0 || var22<=0 || scale<=0 || smooth<0 || fabs(col)>=1) rho=-2;
     break;
     case 124:   /*parsimonious LMC*/
nug11=par[3];
nug22=par[4];
scale11=par[5];
scale22=par[6];
if(nug11<0 || nug22<0 ||  scale11<=0 || scale22<=0   ) rho=-2;
     break;
     case 126:   /*not parsimonious LMC*/
nug11=par[4];
nug22=par[5];
scale11=par[6];
scale22=par[7];
if(nug11<0 || nug22<0 ||  scale11<=0 || scale22<=0   ) rho=-2;
     break;
     case 118:  //bivariate matern with contrainsts
     case 121:
     var11=par[0];
     var22=par[1];
     nug11=par[2];
     nug22=par[3];
     col=par[4];
     scale11=par[5];
     scale22=par[6];
     smoo11=par[7];
     smoo22=par[8];
     scale12=0.5*(scale11+scale22);
     smoo12=0.5*(smoo11+smoo22);
     if(nug11<0 || nug22<0 || var11<=0 || var22<=0 || scale11<=0 || scale12<=0 ||
        scale22<0|| smoo11<0 || smoo12<0 || smoo22<0 || fabs(col)>1) { rho=-2;break; }

       /*
     if(smoo12<(0.5*(smoo11+smoo22)))  {if(col!=0) { rho=-2;break; }} //ok

     if((smoo12==(0.5*(smoo11+smoo22))) && (scale12<=fmin(scale11,scale22)))

     { if(R_pow(col,2)>bi_matern_bounds(scale11,scale22,scale12,smoo11,smoo22,smoo12,0,0)) { rho=-2;break;}} //ok

     if((smoo12==(0.5*(smoo11+smoo22))) && (scale12> fmin(scale11,scale22)) && (scale12<fmax(scale11,scale22)))
     {
      if(R_pow(col,2)>bi_matern_bounds(scale11,scale22,scale12,smoo11,smoo22,smoo12,0,0) &&
         R_pow(col,2)>bi_matern_bounds(scale11,scale22,scale12,smoo11,smoo22,smoo12,-LOW,0) &&
         R_pow(col,2)>bi_matern_bounds(scale11,scale22,scale12,smoo11,smoo22,smoo12,
         sqrt(
         ((2*smoo22+2)*R_pow(scale12*scale11,2)+(2*smoo11+2)*R_pow(scale12*scale22,2)-2*(smoo11+smoo22+2)*R_pow(scale22*scale11,2))/
         ((2*smoo11+2)*R_pow(scale11,2)+(2*smoo22+2)*R_pow(scale22,2)-2*(smoo11+smoo22+2)*R_pow(scale12,2))
         ),0))
         {rho=-2;break;}
     }
     if((smoo12==(0.5*(smoo11+smoo22))) && (scale12>=fmax(scale11,scale22)))
     {
         if(R_pow(col,2)>bi_matern_bounds(scale11,scale22,scale12,smoo11,smoo22,smoo12,0,1))
         {rho=-2;break;}
     }
     if(smoo12>(0.5*(smoo11+smoo22)))  {
        if(R_pow(col,2)>bi_matern_bounds(scale11,scale22,scale12,smoo11,smoo22,smoo12,0,0)&&
        R_pow(col,2)>bi_matern_bounds(scale11,scale22,scale12,smoo11,smoo22,smoo12,-LOW,0))
        {rho=-2;break;}
     }*/
     break;

   case 128:  //bivariate matern
   case 117: //bivariate F_sphere
     var11=par[0];
     var22=par[1];
     nug11=par[2];
     nug22=par[3];
     col=par[4];
     scale11=par[5];
     scale12=par[6];
     scale22=par[7];
     smoo11=par[8];
     smoo12=par[9];
     smoo22=par[10];
     if(nug11<0 || nug22<0 || var11<=0 || var22<=0 || scale11<=0 || scale12<=0 ||
        scale22<=0|| smoo11<=0 || smoo12<=0 || smoo22<=0|| fabs(col)>1) { rho=-2;break; }


     /*
     if((smoo12==(0.5*(smoo11+smoo22))) && (scale12<=fmin(scale11,scale22)))
     { if(R_pow(col,2)>bi_matern_bounds(scale11,scale22,scale12,smoo11,smoo22,smoo12,0,0)) { rho=-2;break;}} //ok

     if((smoo12==(0.5*(smoo11+smoo22))) && (scale12> fmin(scale11,scale22)) && (scale12<fmax(scale11,scale22)))
     {
      if(R_pow(col,2)>bi_matern_bounds(scale11,scale22,scale12,smoo11,smoo22,smoo12,0,0)&&
         R_pow(col,2)>bi_matern_bounds(scale11,scale22,scale12,smoo11,smoo22,smoo12,-LOW,0)&&
         R_pow(col,2)>bi_matern_bounds(scale11,scale22,scale12,smoo11,smoo22,smoo12,
         sqrt(
         ((2*smoo22+2)*R_pow(scale12*scale11,2)+(2*smoo11+2)*R_pow(scale12*scale22,2)-2*(smoo11+smoo22+2)*R_pow(scale22*scale11,2))/
         ((2*smoo11+2)*R_pow(scale11,2)+(2*smoo22+2)*R_pow(scale22,2)-2*(smoo11+smoo22+2)*R_pow(scale12,2))
         ),0))
         {rho=-2;break;}
     }
     if((smoo12==(0.5*(smoo11+smoo22))) && (scale12>=fmax(scale11,scale22)))
     {
         if(R_pow(col,2)>bi_matern_bounds(scale11,scale22,scale12,smoo11,smoo22,smoo12,0,1))
         {rho=-2;break;}
     }

    if(smoo12<(0.5*(smoo11+smoo22)))  {if(col!=0) { rho=-2;break; }} //ok

     if(smoo12>(0.5*(smoo11+smoo22)))  {
         if(R_pow(col,2)>bi_matern_bounds(scale11,scale22,scale12,smoo11,smoo22,smoo12,0,0)&&
           R_pow(col,2)>bi_matern_bounds(scale11,scale22,scale12,smoo11,smoo22,smoo12,-LOW,0))
         {rho=-2;break;}
     }*/
     break;
   case 134:  //bivariate wendhole 1
   //case 69: //bivariate wendhole 2
     var11=par[0];
     var22=par[1];
     nug11=par[2];
     nug22=par[3];
     col=par[4];
     scale11=par[5];
     scale12=par[6];
     scale22=par[7];
     smoo11=par[8];
     smoo12=par[9];
     smoo22=par[10];
    if(nug11<0 || nug22<0 || var11<=0 || var22<=0 || scale11<=0 || scale12<=0 || scale22<=0 || fabs(col)>=1) rho=-2;
     break;
    }
  return rho;
}
// list of spatial and spatial-temporal correlation functions:

double CorFct(int *cormod, double h, double u, double *par, int c11, int c22)
{
  double arg=0.0, col=0.0,R_power=0.0, R_power1=0.0, R_power2=0.0, R_power_s=0.0, R_power_t=0.0, var11=0.0, var22=0.0;
  double rho=0.0, sep=0, scale=0.0, smooth=0.0,smooth_s=0.0,smooth_t=0.0, scale_s=0.0, scale_t=0, x=0, nug11=0.0, nug22=0.0;
  double scale11=0.0, scale22=0.0, scale12=0.0, smoo11=0.0, smoo22=0.0, smoo12=0.0,R_power11=0.0, R_power22=0.0, R_power12=0.0;
  switch(*cormod) // Correlation functions are in alphabetical order
    {
    case 1:// Cauchy correlation function
      R_power1=2;
      R_power2=par[0];
      scale=par[1];
      rho=CorFunCauchy(h, R_power2, scale);
      break;
    case 2:// Matern1
      scale=par[0];
      rho=exp(-h/scale)*(1+h/scale);
      break;
    case 3:// Matern2
      scale=par[0];
       rho=exp(-h/scale)*(1+h/scale+R_pow(h/scale,2)/3);
      break;
    case 4:// Exponential correlation function
      R_power=1;
      scale=par[0];
      rho=CorFunStable(h, R_power, scale);
      break;
    case 5: // Dagum
    R_power1=par[0];
    R_power2=par[1];
    scale=par[2];
    rho=CorFunDagum(h, R_power1, R_power2, scale);
    break;
    case 8: // Generalised Cuachy correlation function
      R_power1=par[0];
      R_power2=par[1];
      scale=par[2];
      rho=CorFunGenCauchy(h, R_power1, R_power2, scale);
      break;
    case 10:// Skarofski correlation function
      scale_s=par[0];
      scale_t=par[1];
      smooth=par[2];
      rho=Shkarofski(h*h, scale_s,scale_t,smooth);
      break;
    case 11://wen0
        R_power=par[0];
        scale=par[1];
    rho=CorFunW0(h,scale,R_power);
        break;
    case 13://wen1
        R_power=par[0];
        scale=par[1];
    rho=CorFunW1(h,scale,R_power);
        break;
    case 12:// Stable correlation function
      R_power=par[0];
      scale=par[1];
      rho=CorFunStable(h, R_power, scale);
      break;
    case 14://  Whittle-Matern correlation function
      scale=par[0];
      smooth=par[1];
      rho=CorFunWitMat(h, scale, smooth);
      break;
    case 15://wen2
        R_power=par[0];
        scale=par[1];
        rho=CorFunW2(h,scale,R_power);
        break;
    case 16: //wave
      scale=par[0];
      rho=CorFunWave(h,scale);
      break;
    case 17://  multiquadric correlation function valid on sphere
        R_power=par[0];
        scale=par[1];
        rho=R_pow(1-R_power/2,2*scale)/R_pow(1+R_pow(R_power/2,2)-R_power*cos(h/REARTH[0]),scale);


    break;
    case 18://  sinsphere correlation function valid on sphere
        R_power=par[0];
        rho=1-R_pow(sin(h/(2*REARTH[0])),R_power);
    break;
/*######*/
    case 21: // hyperg correlation 2 parameters
        R_power=par[0];
        R_power1=par[1];
        scale=par[2];
        smooth=par[3];
  rho=CorFunHyperg2(h,R_power, R_power1, smooth, scale);
        break;
case 22: // hyperg correlation 1 parameter
        R_power=par[0];
        scale=par[1];
        smooth=par[2];
       rho=CorFunHyperg(h,R_power, smooth, scale);
        break;
case 23: // hyperg correlation  parameter with matern my version
        R_power1=1/par[0];
        scale=par[1];
        smooth=par[2];
        sep=exp(  (lgammafn((R_power1+1)/2+smooth)+lgammafn((R_power1+2+1)/2+2*smooth)
                       -lgamma((R_power1+2)/2+smooth)-lgamma(R_power1/2) )/ (1+2*smooth) 
               );
  rho=CorFunHyperg(h,R_power1, smooth, 2*scale*sep);
        break;
 case 30: // hyperg correlation  parameter with matern emery version
        R_power1=1/par[0];
        scale=par[1];
        smooth=par[2];
  rho=CorFunHyperg(h,R_power1, smooth, scale*R_power1);
        break;       
/*######*/
    case 19: // original   Generalised wend correlation
        R_power1=par[0];
        scale=par[1];
        smooth=par[2];
  rho=CorFunW_gen(h, R_power1, smooth, scale);
      break;
    case 26: // original   Generalised wend hole
         R_power=par[0];
        R_power1=par[1];
        scale=par[2];
        smooth=par[3];
  rho=CorFunW_genhole(h, R_power1, smooth, scale, R_power);
        break;
    case 29: // reparametrized original   Generalised wend hole
         R_power=par[0];
        R_power1=par[1];
        scale=par[2];
        smooth=par[3];
          rho=CorFunW_genhole(h, R_power1, smooth-0.5, scale*R_power1, R_power);
        break;
   case 24: // original   kummer
        R_power1=par[0];
        scale=par[1];
        smooth=par[2];
  rho=CorKummer(h, R_power1, smooth, scale);
        break;
    case 25: //    kummer matern
        R_power1=par[0];
        scale=par[1];
        smooth=par[2];
  rho=CorKummer(h, R_power1, smooth, scale*sqrt(2*(R_power1+1)));
        break;

 case 27://  Whittle-Matern correlation function with hole effect
      R_power1=par[0];
      scale=par[1];
      smooth=par[2];
      rho=CorFunWitMathole(h, scale, smooth,R_power1);
      break;
 case 28://  Schoemberg with hole effect;
      scale=par[0];
      rho=Corschoenberg(h, scale);
      break;
    case 6: // Bevilacqua Generalised wend correlation function  "better"  parametrization
        R_power1=1/par[0];
        scale=par[1];
        smooth=par[2];
        sep=exp(  (lgammafn(2*smooth+R_power1+1)-lgammafn(R_power1))/ (1+2*smooth) );
        rho=CorFunW_gen(h, R_power1, smooth,  scale * sep);
        break;  
     case 7: // emery second parametrization
        R_power1=1/par[0];
        scale=par[1];
        smooth=par[2];
        rho=CorFunW_gen(h, R_power1, smooth,  scale * R_power1);
        break;
     case 20://  F_sphere correlation function
      scale=par[0];
      smooth=par[1];
      rho=CorFunSmoke(h, scale, smooth);
      break;

      /*case 20: // Generalised wend correlation function
        R_power1=par[0];
        scale=par[1];
        smooth=par[2];
        rho=(R_pow(1,2*smooth+1)*CorFunW_gen(h, smooth+3+1.5, smooth, 1)-R_pow(0.75,2*smooth+1)*CorFunW_gen(h, smooth+3+1.5, smooth, 0.75) )/
           (R_pow(1,2*smooth+1)-R_pow(0.75,2*smooth+1));
        break;*/
 /***************** spatial tapers****************************/
   /*case 28:// Bohman taper
      rho=CorFunBohman(h,maxdist[0]);
   case 29:// Bohman model
      rho=CorFunBohman(h,par[0]);
        break;*/
  //case 30:// Wendland0 for tap
  //    rho= CorFunW0(h,maxdist[0],2);
  //    break;
   //case 31:// Wendland1 for model
   //   rho=CorFunW0(h,par[0],2);
   /*         break;
    case 32:// Wendland1 for tap
      rho=CorFunW1(h,maxdist[0],3);
      break;
    case 33:// Wendland1 for tap
      rho=CorFunW1(h,par[0],3);
      break;
    case 34:// Wendland2 for tap
      rho=CorFunW2(h,maxdist[0],4);
      break;
    case 35:// Wendland1 for tap
     rho=CorFunW2(h,par[0],4);
    case 38:// phericalfor tap
     rho=CorFunSferical(h, maxdist[0]);
    break;*/
        case 36:// unit taper
        case 37:// unit taper
      rho=1;
      break;
 /***************** end spatial tapers****************************/
      // START non-separable correlation functions:
    case 42:   //Gneiting correlation model as in (14) Gneitint (2002) with tau=1
      R_power_s=par[0];
      R_power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      sep=par[4];//1/(1+exp(-par[4]));
      arg=1+R_pow(u/scale_t, R_power_t);
      rho=exp(-(R_pow(h/scale_s, R_power_s))*R_pow(arg, -0.5*sep*R_power_s))/arg;
      break;
    case 44:// Iaco-Cesare model as in (14) of Gneitint (2006): note that R_power parameters are in [0,1]
      R_power2=par[0];
      R_power_s=par[1];
      R_power_t=par[2];
      scale_s=par[3];
      scale_t=par[4];
      rho=R_pow(1+R_pow(h/scale_s, R_power_s)+R_pow(u/scale_t, R_power_t),-R_power2);
      break;
    case 46://Porcu model as in (4) of Porcu Bevilacqua et al (2010), with beta_1=beta_2=1
      R_power_s=par[0];
      R_power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      sep=par[4];
      if(sep) rho=R_pow(0.5*R_pow(1+R_pow(h/scale_s, R_power_s),sep)+0.5*R_pow(1+R_pow(u/scale_t, R_power_t),sep),-1/sep);
      else rho=R_pow((1+R_pow(h/scale_s, R_power_s))*(1+R_pow(u/scale_t,R_power_t)),-1);
      break;
    case 48:// Stein model as in (16) of JASA (2005) with epsilon=0*/
      R_power_t=par[0];
      scale_s=par[1];
      scale_t=par[2];
      smooth=par[3];
      arg=smooth+R_pow(u, 0.5*R_power_t)/scale_t;
      if(h==0) {rho=1/(R_pow(2, arg)*gammafn(arg+1));}
      else {rho=R_pow(h/scale_s, arg)*bessel_k(h/scale_s, arg, 1)/(R_pow(2, arg)*gammafn(arg+1));}
      break;
    case 50:   //Gneiting correlation with prac ranges "giusto"
      R_power_s=par[0];
      R_power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      sep=par[4];
      arg=1+R_pow(u, R_power_t)/scale_t;
      rho=exp(-(R_pow(h, R_power_s)/scale_s)*R_pow(arg, 0.5*sep*R_power_s))/R_pow(arg,1.5);
      break;
     case 52:// Gneiting correlation model valid on the sphere (equation 8)
      R_power_s=par[0];
      R_power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      sep=par[4];
      arg=1+R_pow(h/scale_s, 2*R_power_s);
      rho=exp(-R_pow(u/scale_t, 2*R_power_t)/(R_pow(arg, sep*R_power_t)))/arg;
      break;
    case 54:// Gneiting correlation model valid on the sphere  (equazione 9)
      R_power_s=par[0];
      R_power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      sep=par[4];
     // arg=1+R_pow(h/scale_s, 2*R_power_s);
     // rho=R_pow(1+R_pow(u/scale_t, 2*R_power_t)/R_pow(arg,R_power_t*sep),-1)/R_pow(arg,1);
            arg=1+R_pow(u/scale_t, 2*R_power_t);
            rho=exp(-R_pow(h/scale_s, R_power_s)*R_pow(arg,R_power_s*sep))/R_pow(arg,1);


      break;
     case 56:    //st sinR_power
        R_power_s=1;
        R_power_t=par[0];
        scale_s=par[1];
        scale_t=par[2];
        arg=R_pow(1+R_pow(u/scale_t,R_power_t),-1);
        sep=cos(h)*arg;
        rho=(exp(R_power_s*sep/scale_s)*(1+R_power_s*sep/scale_s))/((1+R_power_s/scale_s)*exp(R_power_s/scale_s));


      break;
     case 58:  //st multiquaderic
        R_power_s=par[0];
        R_power_t=par[1];
        scale_s=par[2];
        scale_t=par[3];
      arg=R_pow(1+R_pow(u/scale_t,R_power_t),-1);
     rho= R_pow(R_pow(1-R_power_s/2,2)/(1+R_pow(R_power_s/2,2)-R_power_s*arg*cos(h)),scale_s);   // model B2 in the paper  (eq 9 right part)
     break;
    case 60:
      R_power_s=par[0];
      R_power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      sep=par[4];
      arg=R_pow((1+(R_pow(u/scale_t,R_power_t))),0.5);
      if(h<(scale_s/R_pow(arg,sep))) rho=R_pow(1-h*R_pow(arg,sep)/scale_s,1.5+R_power_s)/arg;
      else                   rho=0;
      break;
    case 61:  //no sep gneiting  with temporal matern margin
        R_power_s=par[0];
        R_power=par[1];
        scale_s=par[2];
        scale_t=par[3];
        sep=par[4];
        smooth=par[5];
        arg=1+R_pow(h/scale_s, R_power_s);
        if(u) rho=R_pow(arg,-R_power)*CorFunWitMat(u,scale_t*R_pow(arg,sep/2),smooth);
        else  rho=R_pow(arg,-R_power);

        break;

     case 62:  //no sep gneiting  with spatial matern margin

        R_power_t=par[0];
        R_power=par[1];
        scale_s=par[2];
        scale_t=par[3];
        sep=par[4];
        smooth=par[5];
        arg=1+R_pow(u/scale_t, R_power_t);
        if(h)  rho=R_pow(arg,-R_power)*CorFunWitMat(h,scale_s*R_pow(arg,sep/2),smooth);
        else  rho=R_pow(arg,-R_power);
        break;

        case 87:
          R_power_t=par[0];
        R_power_s=par[1];
        R_power=par[2];
        scale_s=par[3];
        scale_t=par[4];
        sep=par[5];
        smooth=par[6];
        arg=R_pow(1+R_pow(u/scale_t,R_power_t),-1);
        rho=R_pow(arg,R_power)*CorFunW_gen(h,R_power_s,smooth,scale_s*R_pow(arg,sep));
        break;
      case 88:
        R_power_s=par[0];
        R_power=par[1];
        R_power_t=par[2];
        scale_s=par[3];
        scale_t=par[4];
        sep=par[5];
        smooth=par[6];
        arg=R_pow(1+R_pow(h/scale_s,R_power_s),-1);
        rho=R_pow(arg,R_power)*CorFunW_gen(u,R_power_t,smooth,scale_t*R_pow(arg,sep));
        break;
   case 89:
      scale_s=par[0];
      scale_t=par[1];
      smooth_s=par[2];
       smooth_t=par[3];
       sep=par[4];
     rho=Matern_Matern_nosep(h,u,scale_s,scale_t,smooth_s,smooth_t,sep);
       break;

          case 85:
      scale_s=par[0];
      scale_t=par[1];
      R_power_s=par[2];
      R_power_t=par[3];
      smooth_s=par[4];
       smooth_t=par[5];
       sep=par[6];
     rho=GenWend_GenWend_nosep(h,u,scale_s,scale_t, R_power_s,R_power_t,smooth_s,smooth_t,sep);
       break;


    case 63:  //
         R_power_t=par[0];
        R_power_s=par[1];
        R_power=par[2];
        scale_s=par[3];
        scale_t=par[4];
        sep=par[5];
        //arg=R_pow(1+R_pow(u/scale_t,R_power_t/2),-1/(R_power_t/2));
        arg=R_pow(1+R_pow(u/scale_t,R_power_t),-1);
        rho=R_pow(arg,R_power)*CorFunW0(h,scale_s*R_pow(arg,sep),R_power_s);
        break;
          case 64:
        R_power_s=par[0];
        R_power=par[1];
        R_power_t=par[2];
        scale_s=par[3];
        scale_t=par[4];
        sep=par[5];
        //arg=R_pow(1+R_pow(h/scale_s,R_power_s),-1);
        //rho=R_pow(arg,R_power)*CorFunW0(u,scale_t*R_pow(arg,sep),R_power_t);
        arg=R_pow(1+R_pow(h/scale_s,R_power_s),-1);
       /*  arg=R_pow(1+R_pow(h/scale_s,R_power_s),-1);  */
        rho=R_pow(arg,R_power)*CorFunW0(u,scale_t*R_pow(arg,sep),R_power_t);  //2.5+2*0
         //2.5+2*0
        break;
    case 65:
          R_power_t=par[0];
        R_power_s=par[1];
        R_power=par[2];
        scale_s=par[3];
        scale_t=par[4];
        sep=par[5];
        arg=R_pow(1+R_pow(u/scale_t,R_power_t),-1);
        rho=R_pow(arg,R_power)*CorFunW1(h,scale_s*R_pow(arg,sep),R_power_s);
        break;

     case 66:
        R_power_s=par[0];
        R_power=par[1];
        R_power_t=par[2];
        scale_s=par[3];
        scale_t=par[4];
        sep=par[5];
        arg=R_pow(1+R_pow(h/scale_s,R_power_s),-1);
        rho=R_pow(arg,R_power)*CorFunW1(u,scale_t*R_pow(arg,sep),R_power_t); //2.5+2*1
        break;
      case 67:  //
         R_power_t=par[0];
        R_power_s=par[1];
        R_power=par[2];
        scale_s=par[3];
        scale_t=par[4];
        sep=par[5];
        arg=R_pow(1+R_pow(u/scale_t,R_power_t),-1);
        rho=R_pow(arg,R_power)*CorFunW2(h,scale_s*R_pow(arg,sep),R_power_s);
        break;

     case 68:
        R_power_s=par[0];
        R_power_t=par[2];
        R_power=par[1];
        scale_s=par[3];
        scale_t=par[4];
        sep=par[5];
        arg=R_pow(1+R_pow(h/scale_s,R_power_s),-1);
        rho=R_pow(arg,R_power)*CorFunW2(u,scale_t*R_pow(arg,sep),R_power_t); ////2.5+2*2
        break;
            case 69:
          R_power_s=par[0];
          R_power_t=par[1];
          scale_s=par[2];
          scale_t=par[3];

          rho=CorFunW_gen(h, R_power_s, 0, scale_s)*CorFunW_gen(u, R_power_t, 0, scale_t);
          //rho=CorFunW0(h,scale_s,R_power_s)*CorFunW0(u,scale_t,R_power_t);
        break;
        case 70:
              R_power_s=par[0];
          R_power_t=par[1];
          scale_s=par[2];
          scale_t=par[3];
          rho=CorFunW0(h,scale_s,R_power_s)*CorFunW1(u,scale_t,R_power_t);
        break;
        case 71:
              R_power_s=par[0];
          R_power_t=par[1];
          scale_s=par[2];
          scale_t=par[3];
          rho=CorFunW0(h,scale_s,R_power_s)*CorFunW2(u,scale_t,R_power_t);
        break;
        case 72:
              R_power_s=par[0];
          R_power_t=par[1];
          scale_s=par[2];
          scale_t=par[3];
          rho=CorFunW1(h,scale_s,R_power_s)*CorFunW0(u,scale_t,R_power_t);
        break;
        case 73:
          R_power_s=par[0];
          R_power_t=par[1];
          scale_s=par[2];
          scale_t=par[3];
          rho=CorFunW1(h,scale_s,R_power_s)*CorFunW1(u,scale_t,R_power_t);
        break;
          case 74:
                R_power_s=par[0];
          R_power_t=par[1];
          scale_s=par[2];
          scale_t=par[3];
          rho=CorFunW1(h,scale_s,R_power_s)*CorFunW2(u,scale_t,R_power_t);
        break;
          case 75:
                R_power_s=par[0];
          R_power_t=par[1];
          scale_s=par[2];
          scale_t=par[3];
          rho=CorFunW2(h,scale_s,R_power_s)*CorFunW0(u,scale_t,R_power_t);
        break;
          case 76:
                R_power_s=par[0];
          R_power_t=par[1];
          scale_s=par[2];
          scale_t=par[3];
          rho=CorFunW2(h,scale_s,R_power_s)*CorFunW1(u,scale_t,R_power_t);
        break;
        case 77:
              R_power_s=par[0];
          R_power_t=par[1];
          scale_s=par[2];
          scale_t=par[3];
          rho=CorFunW2(h,scale_s,R_power_s)*CorFunW2(u,scale_t,R_power_t);
        break;
          case 78:
            scale_s=par[0];
            scale_t=par[1];
            smooth_s=par[2];
            smooth_t=par[3];
            R_power_s=par[4];
            R_power_t=par[5];
          rho=CorFunW_gen(h, R_power_s, smooth_s, scale_s)*CorFunW_gen(u, R_power_t, smooth_t, scale_t);
        break;
      // END non-separable correlation functions
      // START separable correlation functions:
    case 82:// Exp-Cauchy:
      R_power2=par[0];
      scale_s=par[1];
      scale_t=par[2];
      rho=CorFunStable(h,1,scale_s)*R_pow((1+R_pow(u/scale_t, 2)), -R_power2);
      break;
    case 84:// Double exp:
      scale_s=par[0];
      scale_t=par[1];
      rho=CorFunStable(h,1,scale_s)*CorFunStable(u,1,scale_t);
      break;
     case 96:// prove
        R_power_s=par[0];
          R_power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      rho=exp(-    R_pow(h/scale_s,R_power_s)* exp(-u/scale_t)       );
      break;
    case 86:// Matern Matern
      scale_s=par[0];
      scale_t=par[1];
      smooth_s=par[2];
       smooth_t=par[3];
      rho=CorFunWitMat(h, scale_s,smooth_s)*CorFunWitMat(u, scale_t, smooth_t);
      break;
    case 90:// Matern-Cauchy:
      R_power2=par[0];
      scale_s=par[1];
      scale_t=par[2];
      smooth=par[3];
      if(h==0) {arg=1;}
      else {arg=R_pow(2,1-smooth)/gammafn(smooth)*R_pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth, 1);}
      rho=arg*R_pow((1+R_pow(u/scale_t, 2)),-R_power2);
      break;
    case 92:// Matern-exp:
      scale_s=par[0];
      scale_t=par[1];
      smooth=par[2];
      if(h==0) {arg=1;}
      else {arg=R_pow(2,1-smooth)/gammafn(smooth)*R_pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth, 1);}
      rho=arg*exp(-u/scale_t);
      break;
    case 94:// Stable-stab:
      R_power_s=par[0];
      R_power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      rho=CorFunStable(h,R_power_s,scale_s)*CorFunStable(u,R_power_t,scale_t);
      break;
      /****************************** bivariate ***************************************/

   case 122:       //bivariate sep matern
     var11=par[0];
     var22=par[1];
     nug11=par[2];
     nug22=par[3];
     col=par[4];
     scale=par[5];
     smooth=par[6];
     if((c11==0)&&(c22==0))                   {if(h==0)  {rho=var11+nug11;}
                                              else       {rho=var11*CorFunWitMat(h, scale, smooth);}
                                              break;}
     if((c11==0&&c22==1)||(c11==1&&c22==0))   {
                                          if(h==0) {rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));}
                                          else     {rho=col*sqrt(var11)*sqrt(var22)*CorFunWitMat(h, scale, smooth);}
                                          break;}
     if((c11==1)&&(c22==1))                   {if(h==0)  {rho=var22+nug22;}
                                              else       {rho=var22*CorFunWitMat(h, scale, smooth);}
                                              break;}
        break;

         case 119:       //bivariate F_sphere sep
     var11=par[0];
     var22=par[1];
     nug11=par[2];
     nug22=par[3];
     col=par[4];
     scale=par[5];
     smooth=par[6];
     if((c11==0)&&(c22==0))                   {if(h==0)  {rho=var11+nug11;}
                                              else      { rho=var11*CorFunSmoke(h, scale, smooth);}
                                              break;}
     if((c11==0&&c22==1)||(c11==1&&c22==0))   {
                                          if(h==0) {rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));}
                                          else     {rho=col*sqrt(var11)*sqrt(var22)*CorFunSmoke(h, scale, smooth);}
                                          break;}
     if((c11==1)&&(c22==1))                   {if(h==0)  {rho=var22+nug22;}
                                              else       {rho=var22*CorFunSmoke(h, scale, smooth);}
                                              break;}
        break;


        case 124:   /*parsimonious LMC */
        var11=par[0];
        col=par[1];
        var22=par[2];
        nug11=par[3];
        nug22=par[4];
        scale11=par[5];
        scale22=par[6];
        smoo11=  CorFunStable(h, 1, scale11); //CorFunStable(h, 2, scale11);    //; for environmetrics
        smoo22=CorFunStable(h, 1, scale22); //CorFunWave(h,scale22);        //;      for environmetrics
        if((c11==0)&&(c22==0))                  { if(h==0) {rho=rho+nug11;break;}
                                                  rho=R_pow(var11,2)*smoo11+R_pow(col,2)*smoo22;
                                                 }
        if((c11==0&&c22==1)||(c11==1&&c22==0))  {rho=var11*col*smoo11+var22*col*smoo22;break;}
        if((c11==1)&&(c22==1))                  {if(h==0) {rho=rho+nug22;break;}
                                                 rho=R_pow(col,2)*smoo11+R_pow(var22,2)*smoo22;
                                                 }
        break;

        case 126:   /*not parsimonious LMC */
        var11=par[0];
        col=par[1];
        var22=par[2];
        smooth=par[3];
        nug11=par[4];
        nug22=par[5];
        scale11=par[6];
        scale22=par[7];
        smoo11=CorFunStable(h, 1, scale11);
        smoo22=CorFunStable(h, 1, scale22);
        if((c11==0)&&(c22==0))                  {rho=R_pow(var11,2)*smoo11+R_pow(col,2)*smoo22;
                                                if(h==0) {rho=rho+nug11;break;}}
        if((c11==0&&c22==1)||(c11==1&&c22==0))  {rho=var11*smooth*smoo11+var22*col*smoo22;break;}
        if((c11==1)&&(c22==1))                  {rho=R_pow(smooth,2)*smoo11+R_pow(var22,2)*smoo22;
                                                if(h==0) {rho=rho+nug22;break;}}
        break;
        case 111:       // multi wend(k=1) separable
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        R_power=par[5];
        scale=par[6];
        if((c11==0)&&(c22==0))            {if(h==0)  {rho=var11+nug11;}
        else      { rho=var11*CorFunW0(h,scale,R_power);}break;
                                           }
        if((c11==0&&c22==1)||(c11==1&&c22==0))  {
                                         if(h==0) {rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));}
            else {rho=col*sqrt(var11)*sqrt(var22)*CorFunW0(h,scale,R_power);}break;}
        if((c11==1)&&(c22==1))            {if(h==0)  {rho=var22+nug22;}
        else     { rho=var22*CorFunW0(h,scale,R_power);}
                                           break;}
        break;
        case 113:       // multi wend(k=1) separable
            var11=par[0];
            var22=par[1];
            nug11=par[2];
            nug22=par[3];
            col=par[4];
            R_power=par[5];
            scale=par[6];
            if((c11==0)&&(c22==0))            {if(h==0)  {rho=var11+nug11;}
            else       {rho=var11*CorFunW1(h,scale,R_power);}
                break;}
            if((c11==0&&c22==1)||(c11==1&&c22==0))  {
                                                     if(h==0) {rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));}
                else {rho=col*sqrt(var11)*sqrt(var22)*CorFunW1(h,scale,R_power);}break;}
            if((c11==1)&&(c22==1))            {if(h==0)  {rho=var22+nug22;}
            else      {rho=var22*CorFunW1(h,scale,R_power);}
                break;}
            break;

        case 115:       // multi wend(k=1) separable
            var11=par[0];
            var22=par[1];
            nug11=par[2];
            nug22=par[3];
            col=par[4];
            R_power=par[5];
            scale=par[6];
            if((c11==0)&&(c22==0))            {if(h==0)  {rho=var11+nug11;}
            else       {rho=var11*CorFunW2(h,scale,R_power);}
                break;}
            if((c11==0&&c22==1)||(c11==1&&c22==0))  {
              if(h==0) {rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));}
                                                     else {rho=col*sqrt(var11)*
                                                         sqrt(var22)*CorFunW2(h,scale,R_power);}break;}
            if((c11==1)&&(c22==1))            {if(h==0)  {rho=var22+nug22;}
            else      {rho=var22*CorFunW2(h,scale,R_power);}
                break;}
            break;


        case 112:       // multi wend(k=0)
                var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        R_power11=par[5];
        R_power12=par[6];
        R_power22=par[7];
        scale11=par[8];
        scale22=par[9];
        scale12=par[10];
        if((c11==0)&&(c22==0))   {if(h==0)  {rho=var11+nug11;}
                                  else       {rho=var11*CorFunW0(h,scale11,R_power11);}
                                  break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))  { if(h==0) {rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));}
            else {rho=col*sqrt(var11)*sqrt(var22)*CorFunW0(h,scale12,R_power12);}break;}
        if((c11==1)&&(c22==1))   {if(h==0)  {rho=var22+nug22;}
        else      {rho=var22*CorFunW0(h,scale22,R_power22);}
                                  break;}

        break;
        case 114:       // multi wend(k=1)
            var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        R_power11=par[5];
        R_power12=par[6];
        R_power22=par[7];
        scale11=par[8];
        scale22=par[9];
        scale12=par[10];
        if((c11==0)&&(c22==0))   {if(h==0)  {rho=var11+nug11;}
                                  else       {rho=var11*CorFunW1(h,scale11,R_power11);}
                                  break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))  { if(h==0) {rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));}
            else {rho=col*sqrt(var11)*sqrt(var22)*CorFunW1(h,scale12,R_power12);}break;}
        if((c11==1)&&(c22==1))   {if(h==0)  {rho=var22+nug22;}
        else     { rho=var22*CorFunW1(h,scale22,R_power22);}
                                  break;}

        break;
        case 116:       // multi wend(k=2)
       var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        R_power11=par[5];
        R_power12=par[6];
        R_power22=par[7];
        scale11=par[8];
        scale22=par[9];
        scale12=par[10];
        if((c11==0)&&(c22==0))   {if(h==0)  {rho=var11+nug11;}
                                  else       {rho=var11*CorFunW2(h,scale11,R_power11);}
                                  break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))  { if(h==0) {rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));}
            else {rho=col*sqrt(var11)*sqrt(var22)*CorFunW2(h,scale12,R_power12);}break;}
        if((c11==1)&&(c22==1))   {if(h==0)  {rho=var22+nug22;}
                                  else      {rho=var22*CorFunW2(h,scale22,R_power22);}
                                  break;}
        break;
         case 129:       // multi wend(k=0) contr
          var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        R_power11=par[5];
        R_power12=par[6];
        R_power22=0.5*(R_power11+R_power22);
        scale11=par[7];
        scale22=par[8];
        scale12=0.5*(scale11+scale22);
        if((c11==0)&&(c22==0))   {if(h==0)  {rho=var11+nug11;}
                                  else       {rho=var11*CorFunW0(h,scale11,R_power11);}
                                  break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))  { if(h==0) {rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));}
            else {rho=col*sqrt(var11)*sqrt(var22)*CorFunW0(h,scale12,R_power12);}break;}
        if((c11==1)&&(c22==1))   {if(h==0)  {rho=var22+nug22;}
        else      {rho=var22*CorFunW0(h,scale22,R_power22);}
                                  break;}
        break;
        case 131:       // // multi wend(k=1) contr
            var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        R_power11=par[5];
        R_power22=par[6];
        R_power12=0.5*(R_power11+R_power22);
        scale11=par[7];
        scale22=par[8];
        scale12=0.5*(scale11+scale22);
            if((c11==0)&&(c22==0))   {if(h==0)  rho=var11+nug11;
            else       rho=var11*CorFunW1(h,scale11,R_power11);
                break;}
            if((c11==0&&c22==1)||(c11==1&&c22==0))  { if(h==0) {rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));}
                                                     else {rho=col*sqrt(var11)*
                                                         sqrt(var22)*CorFunW1(h,scale12,R_power12);}break;}
            if((c11==1)&&(c22==1))   {if(h==0)  {rho=var22+nug22;}
            else      {rho=var22*CorFunW1(h,scale22,R_power22);}
                break;}
        break;
        case 120:          // // multi wend(k=2) contr
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        R_power11=par[5];
        R_power22=par[6];
        R_power12=0.5*(R_power11+R_power22);
        scale11=par[7];
        scale22=par[8];
        scale12=0.5*(scale11+scale22);

            if((c11==0)&&(c22==0))   {if(h==0)  {rho=var11+nug11;}
            else      { rho=var11*CorFunW2(h,scale11,R_power11);}
                break;}
            if((c11==0&&c22==1)||(c11==1&&c22==0)) { if(h==0) {rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));}
                                                     else{ rho=col*sqrt(var11)*
                                                         sqrt(var22)*CorFunW2(h,scale12,R_power12);}break;}
            if((c11==1)&&(c22==1))   {if(h==0)  {rho=var22+nug22;}
            else      {rho=var22*CorFunW2(h,scale22,R_power22);}
                break;}
        break;
        case 118:       // full bivariate matern with contraists
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        scale11=par[5];
        scale22=par[6];
        smoo11=par[7];
        smoo22=par[8];
        scale12=0.5*(scale11+scale22);
        smoo12=0.5*(smoo11+smoo22);

        if((c11==0)&&(c22==0))    {if(h==0)  {rho=var11+nug11;}
                                   else      {rho=var11*CorFunWitMat(h, scale11,  smoo11);}
                                   break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0)){ if(h==0) {rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));}
            else {rho=col*sqrt(var11)*sqrt(var22)*CorFunWitMat(h, scale12,  smoo12);}break;}
        if((c11==1)&&(c22==1))   {if(h==0)  {rho=var22+nug22;}
        else      {rho=var22*CorFunWitMat(h, scale22,  smoo22);}
                                  break;}
        break;
        case 121:       // full bivariate F_sphere with contraists
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        scale11=par[5];
        scale22=par[6];
        smoo11=par[7];
        smoo22=par[8];
        scale12=0.5*(scale11+scale22);
        smoo12=0.5*(smoo11+smoo22);
        if((c11==0)&&(c22==0))    {if(h==0)  {rho=var11+nug11;}
        else      {rho=var11*CorFunSmoke(h, scale11,  smoo11);}
                                   break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0)){ if(h==0) {rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));}
            else {rho=col*sqrt(var11)*sqrt(var22)*CorFunSmoke(h, scale12,  smoo12);}break;}
        if((c11==1)&&(c22==1))   {if(h==0)  {rho=var22+nug22;}
        else      {rho=var22*CorFunSmoke(h, scale22,  smoo22);}
                                  break;}
        break;
        case 128:       // full bivariate matern
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        scale11=par[5];
        scale12=par[6];
        scale22=par[7];
        smoo11=par[8];
        smoo12=par[9];
        smoo22=par[10];
        if((c11==0)&&(c22==0))  {if(h==0)  {rho=var11+nug11;}
                                else      {rho=var11*CorFunWitMat(h, scale11,  smoo11);}
                                break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))    {
                                               if(h==0) {rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));}
                                                     else {rho=col*sqrt(var11*
                                                         var22)*CorFunWitMat(h, scale12,  smoo12);}break;}
        if((c11==1)&&(c22==1))  {if(h==0)  {rho=var22+nug22;}
                                 else      {rho=var22*CorFunWitMat(h, scale22,  smoo22);}
                                 break;}
        break;
             case 117:       // full bivariate matern
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        scale11=par[5];
        scale12=par[6];
        scale22=par[7];
        smoo11=par[8];
        smoo12=par[9];
        smoo22=par[10];
        if((c11==0)&&(c22==0))  {if(h==0)  {rho=var11+nug11;}
                                else      {rho=var11*CorFunSmoke(h, scale11,  smoo11);}
                                break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))    {
                                               if(h==0) {rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));}
                                                     else {rho=col*sqrt(var11*
                                                         var22)*CorFunSmoke(h, scale12,  smoo12);}break;}
        if((c11==1)&&(c22==1))  {if(h==0)  {rho=var22+nug22;}
                                 else     { rho=var22*CorFunSmoke(h, scale22,  smoo22);}
                                 break;}
        break;
        /************************************************/
        case 130:       //biv gen wend sep
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        R_power=par[5];
        scale=par[6];
        smooth=par[7];
        if((c11==0)&&(c22==0))            {if(h==0)  {rho=var11+nug11;}
        else       {rho=var11*CorFunW_gen(h, R_power, smooth, scale);}
                                           break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))  { if(h==0){ rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));}
            else {rho=col*sqrt(var11)*sqrt(var22)*CorFunW_gen(h, R_power, smooth, scale);}break;}
        if((c11==1)&&(c22==1))            {if(h==0)  {rho=var22+nug22;}
        else      {rho=var22*CorFunW_gen(h, R_power, smooth, scale);}
                                           break;}
        break;

/************************************************/
              case 132:       //biv gen wend
            var11=par[0];
            var22=par[1];
            nug11=par[2];
            nug22=par[3];
            col=par[4];
            R_power11=par[5];
            R_power12=par[6];
            R_power22=par[7];
            scale11=par[8];
            scale12=par[9];
            scale22=par[10];
            smoo11=par[11];
            smoo12=par[12];
            smoo22=par[13];
        if((c11==0)&&(c22==0))   {if(h==0)  {rho=var11+nug11;}
        else       {rho=var11*CorFunW_gen(h, R_power11, smoo11, scale11);}
                                  break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))  { if(h==0) {rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));}
            else {rho=col*sqrt(var11)*sqrt(var22)*CorFunW_gen(h, R_power12, smoo12, scale12);}break;}
        if((c11==1)&&(c22==1))   {if(h==0)  {rho=var22+nug22;}
        else      {rho=var22*CorFunW_gen(h, R_power22, smoo22, scale22);}
                                  break;}
        break;

/************************************************/

                case 134:       // biv gen wend contr
          var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        R_power11=par[5];
        R_power12=par[6];
        R_power22=0.5*(R_power11+R_power22);
        scale11=par[7];
        scale22=par[8];
        scale12=0.5*(scale11+scale22);
        smoo11=par[9];
        smoo22=par[10];
        smoo12=0.5*(smoo11+smoo22);
        if((c11==0)&&(c22==0))   {if(h==0)  {rho=var11+nug11;}
                                  else       {rho=var11*CorFunW_gen(h, R_power11, smoo11, scale11);}
                                  break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))  { if(h==0) {rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));}
                                                     else {rho=col*sqrt(var11)*
                                                         sqrt(var22)*CorFunW_gen(h, R_power12, smoo12, scale12);}break;}
        if((c11==1)&&(c22==1))   {if(h==0)  {rho=var22+nug22;}
                                  else      {rho=var22*CorFunW_gen(h, R_power22, smoo22, scale22);}
                                  break;}
        break;

    /************************************************/
        case 136:       // matern cauchy
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        scale11=par[5];
        scale12=par[6];
        scale22=par[7];
        smoo11=par[8];
        smoo12=par[9];
        smoo22=par[10];
        R_power22=par[11];


        if((c11==0)&&(c22==0))  {if(h==0)  {rho=var11+nug11;}
                                else      {rho=var11*CorFunWitMat1(h*h, scale11,  smoo11);}
                                break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0)){ if(h==0) {rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));}
                                                else {rho=col*sqrt(var11)*
                                                         sqrt(var22)*Shkarofski(h*h,scale12,scale12,smoo12);}
                                                       //CorFunWitMatCau(h,scale12,smoo12);
                                                       break;}
        if((c11==1)&&(c22==1))  {if(h==0)  {rho=var22+nug22;}
                                 else      {rho=var22*CorFunGenCauchy2(h*h,smoo22,R_power22,scale22);}
                                 break;}
        break;
           case 137:       // gen matern cauchy
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        scale11=par[5];
        scale12=par[6];
        scale22=par[7];
        smoo11=par[8];
        smoo12=par[9];
        smoo22=par[10];
        R_power11=par[11];
        R_power12=par[12];
        R_power22=par[13];

        if((c11==0)&&(c22==0))  {if(h==0)  {rho=var11+nug11;}
        else      {rho=var11*CorFunGenWitMatCau(h, scale11,  smoo11,R_power11);}

                                break;}
if((c11==0&&c22==1)||(c11==1&&c22==0)){ if(h==0) {rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));}
                                                     else {rho=col*sqrt(var11)*sqrt(var22)*CorFunGenWitMatCau(h, scale12,  smoo12,R_power12);}
                                break;}

        if((c11==1)&&(c22==1))  {if(h==0)  {rho=var22+nug22;}
                                 else      {rho=var22*CorFunGenWitMatCau(h, scale22,  smoo22,R_power22);}

                                 break;}
        break;
           case 138:       //  bivariate cauchy
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        scale11=par[5];
        scale12=par[6];
        scale22=par[7];
        R_power11=par[8];
        R_power12=par[9];
        R_power22=par[10];
        smoo11=par[11];
        smoo12=par[12];
        smoo22=par[13];

        if((c11==0)&&(c22==0))  {if(h==0)  {rho=var11+nug11;}
                                else      {rho=var11*CorFunGenCauchy(h, R_power11, smoo11, scale11);}
                                break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))    {
          if(h==0) {rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));}
           else {rho=col*sqrt(var11)*sqrt(var22)*CorFunGenCauchy(h, R_power12, smoo12, scale12);}break;}
        if((c11==1)&&(c22==1))  {if(h==0)  {rho=var22+nug22;}
                                 else      {rho=var22*CorFunGenCauchy(h, R_power22, smoo22, scale22);}
                                 break;}

        break;
        case 139:       //  bivariate stable
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        scale11=par[5];
        scale12=par[6];
        scale22=par[7];
        R_power11=par[8];
        R_power12=par[9];
        R_power22=par[10];
        if((c11==0)&&(c22==0))  {if(h==0)  {rho=var11+nug11;}
                                else      {rho=var11*CorFunStable(h, R_power11, scale11);}
                                break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))    {
          if(h==0) {rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));}
                                                     else {rho=col*sqrt(var11)*sqrt(var22)*CorFunStable(h, R_power12, scale12);}break;}
        if((c11==1)&&(c22==1))  {if(h==0)  {rho=var22+nug22;}
                                 else      {rho=var22*CorFunStable(h, R_power22, scale22);}
                                 break;}

        break;
        /****************************** bivariate  taper ***************************************/
        case 140:       // wendland-gneiting k=0
        if((c11==0)&&(c22==0))                   {rho=CorFunW0(h,dista[0][0],2);break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))   {rho=tapsep[0]*CorFunW0(h,dista[0][1],2);break;}
        if((c11==1)&&(c22==1))                   {rho=CorFunW0(h,dista[1][1],2);break;};
        break;
        case 141:       // wendland-gneiting k=0
        if((c11==0)&&(c22==0))                   {rho=CorFunW0(h,par[0],2);break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))   {rho=par[3]*CorFunW0(h,par[1],2);break;}
        if((c11==1)&&(c22==1))                   {rho=CorFunW0(h,par[2],2);break;};
        break;
        case 142:       // wendland-gneiting k=1
        if((c11==0)&&(c22==0))                    {rho=CorFunW1(h,dista[0][0],3);break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))    {rho=tapsep[0]*CorFunW1(h,dista[0][1],3);break;}
        if((c11==1)&&(c22==1))                    {rho=CorFunW1(h,dista[1][1],3);break;}
        break;
        case 143:       // wendland-gneiting k=1
        if((c11==0)&&(c22==0))                    {rho=CorFunW1(h,par[0],3);break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))    {rho=par[3]*CorFunW1(h,par[1],3);break;}
        if((c11==1)&&(c22==1))                    {rho=CorFunW1(h,par[2],3);break;}
        break;
        case 144:        // wendland-gneiting k=2
        if((c11==0)&&(c22==0))                   {rho=CorFunW2(h,dista[0][0],4);break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))   {rho=tapsep[0]*CorFunW2(h,dista[0][1],4);break;}
        if((c11==1)&&(c22==1))                   {rho=CorFunW2(h,dista[1][1],4);break;}
         break;
        case 145:        // wendland-gneiting k=2
            if((c11==0)&&(c22==0))                   {rho=CorFunW2(h,par[0],4);break;}
            if((c11==0&&c22==1)||(c11==1&&c22==0))   {rho=par[3]*CorFunW2(h,par[1],4);break;}
            if((c11==1)&&(c22==1))                   {rho=CorFunW2(h,par[2],4);break;}
            break;
        case 146:        // Asymmetric taper with wendland-gneiting k=1
        if((c11==0)&&(c22==0)) {  rho=CorFunWend1_tap(h,dista[0][0],0);break;}
        if((c11==0&&c22==1))  {
            rho=1*(tapsep[0]*CorFunWend1_tap(h,dista[0][0],0)+(1-tapsep[0])*CorFunWend1_tap(h,dista[1][1],0));break;}
        if((c11==1&&c22==0))   {
            rho=1*(tapsep[0]*CorFunWend1_tap(h,dista[1][1],0)+(1-tapsep[0])*CorFunWend1_tap(h,dista[0][0],0));break;}
        if((c11==1)&&(c22==1)) {  rho=CorFunWend1_tap(h,dista[1][1],0);break;}
         break;
       case 147:        // Unit

        if((c11==0)&&(c22==0)) {  rho=1;break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))  {rho=1;break;}
        if((c11==1)&&(c22==1)) {  rho=1;break;}
         break;

        /******************************space time taper***************************************/
        case 200:
        rho=CorFunW0(h,maxdist[0],2)*CorFunW0(u,maxtime[0],2);
        break;
        case 201:
        rho=CorFunW0(h,par[0],2)*CorFunW0(u,par[1],2);
        break;
        case 202:
        rho=CorFunW0(h,maxdist[0],2)*CorFunW1(u,maxtime[0],3);
        break;
        case 203:
        rho=CorFunW0(h,par[0],2)*CorFunW1(u,par[1],3);
        break;
        case 204:
        rho=CorFunW0(h,maxdist[0],2)*CorFunW2(u,maxtime[0],4);
        break;
        case 205:
        rho=CorFunW0(h,par[0],2)*CorFunW2(u,par[1],4);
        break;
        case 206:
        rho=CorFunW1(h,maxdist[0],3)*CorFunW0(u,maxtime[0],2);
        break;
        case 207:
        rho=CorFunW1(h,par[0],3)*CorFunW0(u,par[1],2);
        break;
        case 208:
        rho=CorFunW1(h,maxdist[0],3)*CorFunW1(u,maxtime[0],3);
        break;
        case 209:
        rho=CorFunW1(h,par[0],3)*CorFunW1(u,par[1],3);
        break;
        case 210:
        rho=CorFunW1(h,maxdist[0],3)*CorFunW2(u,maxtime[0],4);
        break;
        case 211:
        rho=CorFunW1(h,par[0],3)*CorFunW2(u,par[1],4);
        break;
        case 212:
        rho=CorFunW2(h,maxdist[0],4)*CorFunW0(u,maxtime[0],2);
        break;
        case 213:
        rho=CorFunW2(h,par[0],4)*CorFunW0(u,par[1],2);
        break;
        case 214:
        rho=CorFunW2(h,maxdist[0],4)*CorFunW1(u,maxtime[0],3);
        break;
        case 215:
        rho=CorFunW2(h,par[0],4)*CorFunW1(u,par[1],3);
        break;
        case 216:
        rho=CorFunW2(h,maxdist[0],4)*CorFunW2(u,maxtime[0],4);
        break;
        case 217:
        rho=CorFunW2(h,par[0],4)*CorFunW2(u,par[1],4);
        break;

    /*non separable taper*/
        case 218:  /* non separable temporal adaptive  taper */
        arg=R_pow(1+R_pow(h/maxdist[0],tapsep[1]/2),-tapsep[0]/(tapsep[1]/2));
        rho=R_pow(arg,2.5+2*0)*CorFunW0(u,maxtime[0]*arg,tapsep[2]);
        break;
        case 219:  /* non separable temporal adaptive  taper */
        arg=R_pow(1+R_pow(h/par[0], 1),-par[2]);
        x=u/(arg*par[1]);
        rho=R_pow(arg,3.5)*CorFunW0(x,1,3.5);
        break;
        case 220:  /* non separable temporal adaptive  taper */
        arg=R_pow(1+R_pow(h/maxdist[0],tapsep[1]/2),-tapsep[0]/(tapsep[1]/2));
        rho=R_pow(arg,2.5+2*1)*CorFunW1(u,maxtime[0]*arg,tapsep[2]);
        break;
        case 221:  /* non separable temporal adaptive  taper */
              arg=R_pow(1+R_pow(h/par[0], 1),-par[2]);
        x=u/(arg*par[1]);
        rho=R_pow(arg,5.5)*CorFunW1(x,1,4.5);
        break;
         case 222:  /* non separable temporal adaptive  taper */
        arg=R_pow(1+R_pow(h/maxdist[0],tapsep[1]/2),-tapsep[0]/(tapsep[1]/2));
        rho=R_pow(arg,2.5+2*2)*CorFunW2(u,maxtime[0]*arg,tapsep[2]);
        break;
           case 223:  /* non separable temporal adaptive  taper */
              arg=R_pow(1+R_pow(h/par[0], 1),-par[2]);
        x=u/(arg*par[1]);
        rho=R_pow(arg,7.5)*CorFunW2(x,1,5.5);
        break;


        /*non separable taper*/
        case 224:  /* non separable temporal adaptive  taper */

           R_power_t=par[0];
        R_power=par[1];
        scale_s=par[2];
        scale_t=par[3];
        sep=par[4];
        arg=R_pow(1+R_pow(u/scale_t,R_power_t/2),-sep/(R_power_t/2));
        rho=R_pow(arg,2.5+2*0)*CorFunW0(h,scale_s*arg,R_power);
        break;
         case 225:  /* non separable temporal adaptive  taper */
        arg=R_pow(1+R_pow(u/par[1], 1),-par[2]);
        x=h/(arg*par[0]);
        rho=R_pow(arg,3.5)*CorFunW0(x,1,3.5);
        break;
        case 226:  /* non separable temporal adaptive  taper */
        arg=R_pow(1+R_pow(u/scale_t,tapsep[1]),-tapsep[0]/(tapsep[1]/2));
        rho=R_pow(arg,2.5+2*1)*CorFunW0(h/maxdist[0],arg,tapsep[2]);
        break;
        case 227:  /* non separable temporal adaptive  taper */
              arg=R_pow(1+R_pow(u/par[1], 1),-par[2]);
        x=h/(arg*par[0]);
        rho=R_pow(arg,5.5)*CorFunW1(x,1,4.5);
        break;
         case 228:  /* non separable temporal adaptive  taper */
        arg=R_pow(1+R_pow(u/scale_t,tapsep[1]),-tapsep[0]/(tapsep[1]/2));
        rho=R_pow(arg,2.5+2*2)*CorFunW0(h/maxdist[0],arg,tapsep[2]);
        break;
           case 229:  /* non separable temporal adaptive  taper */
              arg=R_pow(1+R_pow(u/par[1], 1),-par[2]);
        x=h/(arg*par[1]);
        rho=R_pow(arg,7.5)*CorFunW2(x,1,5.5);
        break;
         case 230:  //unit  space time taper
            rho=1;
            break;
        /******************************end space time taper***************************************/
        // END separable correlation functions:

    }
  return rho;
}
// Cauhcy class of correlation models:
double CorFunCauchy(double lag, double R_power2, double scale)
{
  double rho=0.0;
  // Computes the correlation:
  rho=R_pow((1+R_pow(lag/scale,2)),-R_power2/2);
  return rho;
}
// Dagum:
double CorFunDagum(double lag, double R_power1, double R_power2, double scale)
{
    double rho=0.0;
    rho=1-R_pow(R_pow(lag/scale,R_power1)/(1+R_pow(lag/scale,R_power1)), R_power2/R_power1);
    return rho;
}
// Generalised Cauhcy class of correlation models:
double CorFunGenCauchy(double lag, double R_power1, double R_power2, double scale)
{
  double rho=0.0;
  rho=R_pow((1+R_pow(lag/scale,R_power1)), -R_power2/R_power1);
  return rho;
}

// Generalised Cauhcy class of correlation models:
double CorFunGenCauchy2(double lag, double R_power1, double R_power2, double scale)
{
  double rho=0.0;
  rho=R_pow((1+R_pow(lag,R_power2)/scale), -R_power1);
  return rho;
}


double CorFunWitMatCau(double h, double scale12,double smo12)
{
double C,B,dd=h*h;
B=1/bessel_k(1,smo12,1);
C=R_pow(1+dd/(scale12),-smo12/2)*bessel_k(sqrt(1+dd/(scale12)),smo12,1);
return(C*B);
}



double CorFunGenWitMatCau(double h, double scale,double smoo,double beta)
{
double C,B,dd=h*h;
B=1/bessel_k(sqrt(beta/scale),smoo,1);
C=R_pow(1+dd/beta,-smoo/2)*bessel_k(sqrt((dd+beta)/scale),smoo,1);
return(C*B);
}



// Stable class of correlation models:
double CorFunStable(double lag, double R_power, double scale)
{
  double rho=0.0;
  // Computes the correlation:
  rho=exp(-R_pow(lag/scale,R_power));
  return rho;
}
// Double Stable class of correlation models:
/*double CorFunDobStable(double lag, double R_power_s, double R_power_t, double scale_s, double scale_t, double tsep)
{
  double rho=0.0;
  // Computes the correlation:
  rho=exp(-R_pow(lag/scale_s,R_power_s)-R_pow(tsep/scale_t,R_power_t));
  return rho;
  }*/
// Sferical class of correlation models:
double CorFunSferical(double lag, double scale)
{
    double rho=0.0,x=0;
    x=lag/scale;
   if(x<=1) {rho=  R_pow(1-x,2)*(1+x/2);}   else {rho=0;}
   // if(x<=1) {rho=1-1.5 * x+ 0.5 * R_pow(x,3);}  else {rho=0;}
  return rho;
}

// Wave  correlation model:
double CorFunWave(double lag, double scale)
{
  double rho=0.0;
  if(lag==0) { rho=1;}
  else       { rho=(scale/lag)*sin(lag/(scale));}
  return rho;
}

// F_sphere class of correlation models:
double CorFunSmoke(double lag, double scale, double smooth)
{
    if (lag == 0.0) return 1.0;
    double iscale = 1.0 / scale;
    double a = 0.5 + smooth;
    double coslag = cos(lag);
    double one_minus_coslag = 1.0 - coslag;

    // Calcolo log-gamma solo una volta per ogni termine
    double logg1 = lgammafn(iscale + a);
    double logg2 = lgammafn(iscale + smooth);
    double logg3 = lgammafn(2.0 / scale + a);
    double logg4 = lgammafn(smooth);
    double pow_term = R_pow(one_minus_coslag, smooth);
    double hyperg = hypergeo(1.0 / scale + a, 1.0 / scale + smooth, 2.0 / scale + a, coslag);
    double log_coeff = logg1 + logg2 - logg3 - logg4;
    double rho = exp(log_coeff) * pow_term * hyperg;
    return rho;
}


// Whittle=matern class of correlation models:
double CorFunWitMat(double lag, double scale, double smooth)
{
  
    if (lag <= 0.0) return 1.0;
      double a = lag / scale;
    if (smooth == 0.5)  return exp(-a);
    if (smooth == 1.5)  return exp(-a) * (1.0 + a);
    if (smooth == 2.5)  return exp(-a) * (1.0 + a + (a * a) / 3.0);
    if (smooth == 3.5)  return exp(-a) * (1.0 + a + 0.4 * a * a + (a * a * a) / 15.0);

    // Caso generico: formula logaritmica stabile
    double log_besselk = log(bessel_k(a, smooth, 2)); // K_nu(a), con scaling esponenziale interno
    double log_num = smooth * log(a) + log_besselk - a;
    double log_den = (smooth - 1.0) * log(2.0) + lgammafn(smooth);

    return exp(log_num - log_den);
}



double Matern_Matern_nosep(double h, double u,
                           double scale_s, double scale_t,
                           double smooth_s, double smooth_t,
                           double sep)
{

    double C_s = CorFunWitMat(h, scale_s, smooth_s);
    double C_t = CorFunWitMat(u, scale_t, smooth_t);

    if (sep <= 0.0)            /* caso separabile */
        return C_s * C_t;

    /* fattore di interazione non-separabile */
    double denom = 1.0 + sep * h * h * u * u/(1+sep);
    return C_s * C_t / denom;
}




double GenWend_GenWend_nosep(double h, double u,
                           double scale_s, double scale_t,
                            double power_s, double power_t,
                           double smooth_s, double smooth_t,
                           double sep)
{

    double C_s = CorFunW_gen(h, scale_s, power_s, smooth_s);
    double C_t = CorFunW_gen(u, scale_t, power_t, smooth_t);

    if (sep <= 0.0)            /* caso separabile */
        return C_s * C_t;

    /* fattore di interazione non-separabile */
    double denom = 1.0 + sep * h * h * u * u/(1+sep);
    return C_s * C_t / denom;
}




double Shkarofski(double lag, double a,double b, double k)
{
double corr=0.0;
if(a==0 && k>0) return( R_pow(1+sqrt(lag/b),-2*k));
if(b==0 && k<0) return( R_pow(2,1+k) * R_pow(gammafn(-k),-1)  *
                           R_pow(sqrt(lag/a),-k) * bessel_k(sqrt(lag/a),k,1));

corr=R_pow(1+lag/b,-k/2)*bessel_k(sqrt((b+lag)/a),k,1)/bessel_k(sqrt(b/a),k,1);
return(corr);
}

double CorFunWitMat1(double lag, double scale, double smooth)
{
  double rho=0.0;
  double bb=sqrt(lag/scale);
  // Computes the correlation:
  if(lag<=0) rho=1;
  else  rho=(R_pow(2,smooth+1)*R_pow(bb,-smooth)*bessel_k(bb,smooth,1))/(gammafn(-smooth));
  return rho;
}
double CorFunBohman(double lag,double scale)
{
  double rho=0.0,x=0;
  x=lag/scale;
  if(x<=1) {
       if (x>0) rho=(1-x)*(sin(2*M_PI*x)/(2*M_PI*x))+(1-cos(2*M_PI*x))/(2*M_PI*M_PI*x);
       else   rho=1;}
  else rho=0;
  return rho;
}

/* wendland function alpha=0*/
double CorFunW0(double lag, double scale, double smoo)
{

    double x = lag / scale;
     if (x <= 0.0) return 1;  
    if (x >= 1.0) return 0.0;    
    return R_pow(1.0 - x, smoo);
}
/* wendland function alpha=1*/
double CorFunW1(double lag, double scale, double smoo)
{
     double x = lag / scale;
     if (x <= 0.0) return 1;  
     if (x >= 1.0) return 0.0; 
     double s1 = smoo + 1.0;
     return R_pow(1.0 - x, s1) * (1.0 + s1 * x);    
}

/* wendland function alpha=2*/
double CorFunW2(double lag, double scale, double smoo)
{
    double x = lag / scale;
     if (x <= 0.0) return 1;  
     if (x >= 1.0) return 0.0; 
        double s2 = smoo + 2.0;
        double x2 = x * x;
        double s = smoo;
        double poly = 3.0 + x * (3.0 * s + 6.0) + x2 * (s * s + 4.0 * s + 3.0);
        return R_pow(1.0 - x, s2) * poly / 3.0;
}


double CorFunHyperg2(double lag, double R_power, double R_power1, double smooth, double scale)
{
    const double d = 2.0;
    const double epsilon = 1e-14;

    double x = lag / scale;
    if (x <= epsilon) return 1.0;
    if (x >= 1.0)     return 0.0;

    double beta = R_power;
    double gamma = R_power1;
    double a = beta - smooth;
    double b = gamma - smooth;
    double c = beta - smooth + gamma - d / 2.0;

    double logg1 = lgammafn(beta - d / 2.0);
    double logg2 = lgammafn(gamma - d / 2.0);
    double logg3 = lgammafn(c);
    double logg4 = lgammafn(smooth - d / 2.0);

    double one_minus_x2 = 1.0 - x * x;
    double pow_term = R_pow(one_minus_x2, c - 1.0);  // usa R_pow se preferisci la tua funzione

    double hyperg = hypergeo2(a, b, c, one_minus_x2);

    double rho = exp(logg1 + logg2 - logg3 - logg4) * pow_term * hyperg;
    return rho;
}

double CorFunHyperg(double lag, double R_power, double smooth, double scale)
{
    const double d = 2.0;
    const double epsilon = 1e-14;

    double x = lag / scale;
    if (x <= epsilon) return 1.0;
    if (x >= 1.0)     return 0.0;

    double x2 = x * x;
    double one_minus_x2 = 1.0 - x2;

    double rp_half = R_power / 2.0;
    double rp_d_half_smooth = (R_power + d) / 2.0 + smooth;
    double rp_d1_half_2smooth = R_power + (d + 1.0) / 2.0 + 2.0 * smooth;

    // log-gamma
    double logg1 = lgammafn(smooth + (R_power + 1.0) / 2.0);
    double logg2 = lgammafn(2.0 * smooth + (d + R_power + 1.0) / 2.0);
    double logg3 = lgammafn(rp_d1_half_2smooth);
    double logg4 = lgammafn(smooth + 0.5);

    double exponent = R_power + (d - 1.0) / 2.0 + 2.0 * smooth;
    double log_term1 = log(one_minus_x2) * exponent;

    double hyperg = hypergeo2(rp_half, rp_d_half_smooth, rp_d1_half_2smooth, one_minus_x2);

    double log_rho = logg1 + logg2 - logg3 - logg4 + log_term1;
    return exp(log_rho) * hyperg;
}


/* Optimized Wendland correlation function */
/* generalized wendland function*/


double CorFunW_gen(double lag, double R_power1, double smooth, double scale)
{

    double x = lag / scale;

    if (fabs(x) < DBL_EPSILON)
        return 1.0;

    if (x >= 1.0)
        return 0.0;

    if (smooth == 0) {
        return R_pow(1.0 - x, R_power1);
    }

    if (smooth == 1) {
        return R_pow(1.0 - x, R_power1 + 1.0) * (1.0 + x * (R_power1 + 1.0));
    }

    if (smooth == 2) {
        return R_pow(1.0 - x, R_power1 + 2.0) *
               (1.0 + x * (R_power1 + 2.0) + 
               x * x * (R_power1 * R_power1 + 4.0 * R_power1 + 3.0) / 3.0);
    }

    if (smooth == 3) {
        double x2 = x * x;
        double x3 = x2 * x;
        return R_pow(1.0 - x, R_power1 + 3.0) *
               (1.0 +
                x * (R_power1 + 3.0) +
                x2 * (2.0 * R_power1 * R_power1 + 12.0 * R_power1 + 15.0) / 5.0 +
                x3 * (R_power1 * R_power1 * R_power1 + 9.0 * R_power1 * R_power1 + 23.0 * R_power1 + 15.0) / 15.0);
    }

    // General case using hypergeometric function
    return exp(lgammafn(smooth) + lgammafn(2.0 * smooth + R_power1 + 1.0) -
               (lgammafn(2.0 * smooth) + lgammafn(smooth + R_power1 + 1.0))) *
           R_pow(2.0, -R_power1 - 1.0) *
           R_pow(1.0 - x * x, smooth + R_power1) *
           hypergeo2(R_power1 / 2.0, (R_power1 + 1.0) / 2.0,
                     smooth + R_power1 + 1.0, 1.0 - x * x);
}


/* kummer function*/
double CorKummer(double lag,double R_power,double smooth,double scale)  // mu alpha beta
{
  double rho=0.0,x=0.0;
    x=lag/scale;
    if(x<1e-32) rho=1;
    else
    rho=exp((lgammafn(smooth+R_power)-(lgammafn(smooth))))*kummer(R_power,1-smooth,0.5*x*x);
    return(rho);
}    



double CorFunWitMathole(double lag, double scale, double smooth, double R_power1)
{
    double rho = 0.0;
    double d = 2.0;
    double x = lag / scale;
    
    // Early return for zero or very small lag
    if (x < 1e-32) {
        return 1.0;
    }
    
    int k = (int)R_power1;
    
    // Base case: call the simpler function if k is 0
    if (k == 0) {
        return CorFunWitMat(lag, scale, smooth);
    }
    
    // Special case for nu = 0.5 + integer
    if (fabs(smooth - 0.5 - floor(smooth - 0.5)) < 1e-10) {
        int nu_int = (int)(smooth - 0.5);
        
        for (int q = 0; q <= k; q++) {
            for (int r = 0; r <= fmax(0, q - 1); r++) {
                double delta = 0.0;
                
                for (int s = 0; s <= q - r; s++) {
                    for (int t = 0; t <= nu_int; t++) {
                        // Calculate pochammer symbols carefully
                        double poch1 = 1.0;  // poch(smooth + 0.5 - t, 2 * t)
                        for (int i = 0; i < 2 * t; i++) {
                            poch1 *= (smooth + 0.5 - t + i);
                        }
                        
                        double poch2 = 1.0;  // poch(smooth + 0.5 - t - s, s)
                        for (int i = 0; i < s; i++) {
                            poch2 *= (smooth + 0.5 - t - s + i);
                        }
                        
                        // Calculate combination term
                        double comb = 1.0;
                        for (int i = 1; i <= q - r; i++) {
                            comb *= i;
                        }
                        for (int i = 1; i <= s; i++) {
                            comb /= i;
                        }
                        for (int i = 1; i <= (q - r - s); i++) {
                            comb /= i;
                        }
                        
                        // Power of -1
                        double sign = ((q - r - s) % 2 == 0) ? 1.0 : -1.0;
                        
                        // Calculate term
                        double term = sqrt(M_PI) * poch1 * poch2 * sign * comb *
                                     exp(-(smooth - 0.5 + t) * log(2) + 
                                         (smooth - 0.5 - s - t) * log(x) - 
                                         lgammafn(smooth) - 
                                         lgammafn(t + 1) - 
                                         (q - r) * log(scale));
                                         
                        delta += term;
                    }
                }
                
                // Calculate outer pochammer symbols
                double poch_k_q = 1.0;  // poch(k - q + 1, q)
                for (int i = 0; i < q; i++) {
                    poch_k_q *= (k - q + 1 + i);
                }
                
                double poch_q_r = 1.0;  // poch(q, r)
                for (int i = 0; i < r; i++) {
                    poch_q_r *= (q + i);
                }
                
                double poch_q_r2 = 1.0;  // poch(q - r, r)
                for (int i = 0; i < r; i++) {
                    poch_q_r2 *= (q - r + i);
                }
                
                double poch_d2_q = 1.0;  // poch(d/2, q)
                for (int i = 0; i < q; i++) {
                    poch_d2_q *= (d/2 + i);
                }
                
                // Sign for r
                double sign_r = (r % 2 == 0) ? 1.0 : -1.0;
                
                // Calculate q, r term
                rho += sign_r * poch_k_q * poch_q_r * poch_q_r2 / 
                      (R_pow(2.0, q + r) * gammafn(q + 1) * gammafn(r + 1) * poch_d2_q) * 
                      R_pow(lag, q - r) * exp(-x) * delta;
            }
        }
    } 
    else { 
        // General case
        for (int q = 0; q <= k; q++) {
            for (int r = 0; r <= fmax(0, q - 1); r++) {
                double delta = 0.0;
                
                for (int s = 0; s <= q - r; s++) {
                    // Calculate pochammer for s
                    double poch_smooth_s = 1.0;  // poch(smooth + 1 - s, s)
                    for (int i = 0; i < s; i++) {
                        poch_smooth_s *= (smooth + 1 - s + i);
                    }
                    
                    for (int t = 0; t <= q - r - s; t++) {
                        // Calculate combination term
                        double comb = 1.0;
                        for (int i = 1; i <= q - r; i++) {
                            comb *= i;
                        }
                        for (int i = 1; i <= s; i++) {
                            comb /= i;
                        }
                        for (int i = 1; i <= t; i++) {
                            comb /= i;
                        }
                        for (int i = 1; i <= (q - r - s - t); i++) {
                            comb /= i;
                        }
                        
                        // Calculate term
                        delta += comb * poch_smooth_s * 
                                R_pow(-0.5, q - r - s) * 
                                R_pow(x, smooth - s) * 
                                bessel_k(x, smooth + 2 * t + r + s - q, 1);
                    }
                }
                
                // Calculate pochammer symbols for outer terms
                double poch_k_q = 1.0;  // poch(k - q + 1, q)
                for (int i = 0; i < q; i++) {
                    poch_k_q *= (k - q + 1 + i);
                }
                
                double poch_q_r = 1.0;  // poch(q, r)
                for (int i = 0; i < r; i++) {
                    poch_q_r *= (q + i);
                }
                
                double poch_q_r2 = 1.0;  // poch(q - r, r)
                for (int i = 0; i < r; i++) {
                    poch_q_r2 *= (q - r + i);
                }
                
                double poch_d2_q = 1.0;  // poch(d/2, q)
                for (int i = 0; i < q; i++) {
                    poch_d2_q *= (d/2 + i);
                }
                
                // Sign for r
                double sign_r = (r % 2 == 0) ? 1.0 : -1.0;
                
                // Calculate q, r term
                rho += R_pow(x, q - r) * sign_r * poch_k_q * poch_q_r * poch_q_r2 /
                      (R_pow(2.0, q + r) * gammafn(q + 1) * gammafn(r + 1) * poch_d2_q) * delta;
            }
        }
        
        // Apply final constant factor
        rho = rho * R_pow(2.0, 1 - smooth) / gammafn(smooth);
    }
    
    return rho;
}




/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
double CorFunW_genhole(double lag, double R_power1, double smooth, double scale, double kk)
{
    const double d = 2.0;
    double x = lag / scale;
    if (x < 1e-32) return 1.0;

    double mu = R_power1;
    int k = (int)kk;
    if (k == 0) {
        return CorFunW_gen(lag, mu, smooth, scale);
    } else if (x <= 1.0) {
        double alpha = smooth + (d + 1.0) / 2.0 + k;
        double beta  = alpha + mu / 2.0;
        double gama  = beta + 0.5;
        double x2 = x * x;
        int n;
        double cc = 1.0 + d / 2.0 + k;
        double A = 0.0;
        for (n = 0; n <= k; n++) {
            double sign = (n % 2 == 0) ? 1.0 : -1.0;
            double poch1 = poch(cc - beta, n);
            double poch2 = poch(cc - gama, n);
            double poch3 = poch(cc - alpha, n);
            double poch4 = poch(-n, n);
            double a1 = sign * poch1 * poch2 / (poch4 * poch3);
            double a2 = gammafn(n + 1.0) * gammafn(k - n + 1.0) / gammafn(k + 1.0);
            double a3 = R_pow(x, 2 * n);
            double a4 = hypergeo2(cc - beta + n, cc - gama + n, cc - alpha + n, x2);
            A += (a1 * a3 * a4) / a2;
        }
        double uu = d / 2.0 + k;
        double B1 = gammafn(beta - uu) * gammafn(gama - uu) * gammafn(d / 2.0) * gammafn(uu - alpha);
        double B2 = gammafn(uu) * gammafn(alpha - uu) * gammafn(beta - alpha) * gammafn(gama - alpha);
        double B = B1 / B2;
        double tt = 1.0 + alpha;
        double C = 0.0;
        for (n = 0; n <= k; n++) {
            double sign = ((n + k) % 2 == 0) ? 1.0 : -1.0;
            double c1 = sign * poch(1.0 - alpha, k - n) * poch(1.0 + alpha - beta, n) * poch(1.0 + alpha - gama, n) / poch(1.0 + alpha - uu, n);
            double c2 = gammafn(n + 1.0) * gammafn(k - n + 1.0) / gammafn(k + 1.0);
            double c3 = R_pow(x, 2.0 * alpha - d - 2.0 * (k - n));
            double c4 = hypergeo2(tt - beta + n, tt - gama + n, tt - d / 2.0 - k + n, x2);
            C += (c1 * c3 * c4) / c2;
        }
        return A + B * C;
    } else {
        return 0.0;
    }
}
/*******************************************************************************/
double Corschoenberg(double lag,double scale)
{
double rho=0.0,x=0.0;
double d=2;
    x=lag/scale;
    if(x<1e-32) {rho=1; return(rho);}
double d2=d*0.5;
rho=gammafn(d2)*R_pow(x*0.5,1-d2)*bessel_j(x,d2-1);
return(rho);
}    





double CorFunWend0_tap(double lag,double scale,double smoo)
{
    double rho=0.0,x=0;
    x=lag/scale;
    if(x<=1) rho=  R_pow(1-x,smoo+3);
    else rho=0;
    return rho;
}
double CorFunWend1_tap(double lag,double scale,double smoo)
{
    double rho=0.0,x=0;
    x=lag/scale;
    if(x<=1) rho=R_pow(1-x,smoo+5)*(1+(smoo+5)*x);
    else rho=0;
    return rho;
}
double CorFunWend2_tap(double lag,double scale,double smoo)
{
    double rho=0.0,x=0;
    x=lag/scale;
    if(x<=1) rho=R_pow(1-x,smoo+7)*(1+(smoo+7)*x+R_pow(smoo+7,2)*R_pow(x,2)/3);
    else rho=0;
    return rho;
}
double CorFunWend1(double lag,double scale)
{
  double rho=0.0,x=0;
    x=lag/scale;
  if(x<=1) rho=R_pow(1-x,2)*(1+0.5*x);
  else rho=0;
  return rho;
}
double CorFunWend2(double lag,double scale)
{
  double rho=0.0,x=0;
     x=lag/scale;
  if(x<=1) rho=R_pow(1-x,4)*(1+4*x);
  else rho=0;
  return rho;
}
double CorFunWend3(double lag,double scale)
{
  double rho=0.0,x=0;
     x=lag/scale;
  if(x<=1) rho=R_pow(1-x,6)*(1+6*x+(35/3)*x*x);
  else rho=0;
  return rho;
}
double CorFunWend5(double lag,double scale)
{
    double rho=0.0,x=0;
    x=lag/scale;
    if(x<=1) rho=R_pow(1-x,7)*(1+7*x+(48/3)*x*x);
    else rho=0;
    return rho;
}

double CorFunWend4(double lag,double scale)
{
    double rho=0.0,x=0;
    x=lag/scale;
    if(x<=1) rho=R_pow(1-x,5)*(1+5*x);
    else rho=0;
    return rho;
}
double CorFunWendhole1(double lag,double scale)
{
  double rho=0.0,x=0;
     x=lag/scale;
  if(x<=1) rho=R_pow(1-x,4)*(1-3*x);
  else rho=0;
  return rho;
}
double CorFunWendhole2(double lag,double scale)
{
  double rho=0.0,x=0;
     x=lag/scale;
  if(x<=1) rho=R_pow(1-x,5)*(1+4*x-18*R_pow(x,2));
  else rho=0;
  return rho;
}
double CorFunWendhole3(double lag,double scale)
{
  double rho=0.0,x=0;
     x=lag/scale;
  if(x<=1) rho=R_pow(1-x,6)*(1+5*x-3.666667*R_pow(x,2)-70.36111*R_pow(x,3));
  else rho=0;
  return rho;
}
double CorFunWendhole(double lag,double scale)
{
  double rho=0.0,x=0;
     x=lag/scale;
  if(x<=1) rho=R_pow(1-x,5)*(1+5*x-27*R_pow(x,2));
  else rho=0;
  return rho;
}
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/****************** SPATIAL CORRELATION MATRIX (upper trinagular) *******************************/
/************************************************************************************************/
/************************************************************************************************/
// Computation of the upper (lower) triangular spatial correlation matrix: spatial case
void CorrelationMat2(double *rho,double *coordx, double *coordy,double *coordz, double *coordt,  int *cormod,
 double *nuis, double *par,double *radius,int *ns, int *NS)
{
  int i=0,j=0,h=0;// check the paramaters range:
  double dd=0.0;
     for(i=0;i<(ncoord[0]-1);i++){
        for(j=(i+1);j<ncoord[0];j++){
        dd=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],coordz[i],coordz[j],*REARTH);
    rho[h]=CorFct(cormod,dd,0,par,0,0);

       h++;
    }}
  return;
}
// Computation of the upper (lower) triangular spatial discrete  models
void CorrelationMat_dis2(double *rho,double *coordx, double *coordy,double *coordz, double *coordt,  int *cormod, double *mean,
        int *nn,double *nuis, double *par,double *radius, int *ns, int *NS,int *model)
{
    int i=0,j=0,h=0;// check the paramaters range:
    double psj=0.0,dd=0.0,ai=0.0,aj=0.0,p1=0.0,p2=0.0,p=0,corr=0.0,p00=0,p11=0,mui=0.0,muj=0.0;
    for(i=0;i<(ncoord[0]-1);i++){
      for(j=(i+1);j<ncoord[0];j++){

        dd=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],coordz[i],coordz[j],*REARTH);
        corr=CorFct(cormod,dd,0,par,0,0);

   if(*model==14||*model==16||*model==2||*model==11||*model==45){

     ai=mean[i];aj=mean[j];

     psj=pbnorm22(ai,aj,(1-nuis[0])*corr);
     p1=pnorm(ai,0,1,1,0); p2=pnorm(aj,0,1,1,0);
      if(*model==2||*model==11)    
      {
         rho[h]=fmin2(nn[i],nn[j])*(psj-p1*p2);} 
      if(*model==14)       rho[h]=(psj-p1*p2)/((-psj+p1+p2)*p1*p2);
      if(*model==16)       rho[h]=cov_binom_neg(nn[0],psj,p1,p2);
      if(*model==45)     //BinomialNegZINB
      {
           p=pnorm(nuis[2],0,1,1,0);
           p00=pbnorm22(nuis[2],nuis[2],(1-nuis[1])*corr);
           p11=1-2*p+p00;
           dd=cov_binom_neg(nn[0],psj,p1,p2);
           rho[h]=p11*dd +  (R_pow(nn[0],2)*(1-p1)*(1-p2)/(p1*p2)) *(p11-R_pow((1-p),2));
          /// rho[h]=p11*dd +  (R_pow(nn[0],2)* p1*p2 /( (1-p1)*(1-p2) )) *(p11-R_pow((1-p),2));
         }
      }
     if(*model==30||*model==36) // Poisson
       {
           ai=exp(mean[i]);aj=exp(mean[j]);
           rho[h]=sqrt(ai*aj)*corr_pois((1-nuis[0])*corr,ai, aj); // it's the  covariance
       }
 if(*model==46||*model==47) // Poisson gamma
       {
           mui=exp(mean[i]);
           muj=exp(mean[j]);
           ai=mui*(1+mui/nuis[1]);
           aj=muj*(1+muj/nuis[1]);
           rho[h]=sqrt(ai*aj)*corr_pois_gen((1-nuis[0])*corr,mui, muj, nuis[1]); // it's the  covariance
       }
  if(*model==43||*model==44) // poisson inflado
       {
           ai=exp(mean[i]);aj=exp(mean[j]);
           dd=sqrt(ai*aj)*corr_pois((1-nuis[0])*corr,ai, aj); // it's the  covariance poisson
           p=pnorm(nuis[2],0,1,1,0);
           psj=pbnorm22(nuis[2],nuis[2],(1-nuis[1])*corr);
           p1=1-2*p+psj;
          rho[h]=p1*dd +  ai*aj*(p1-R_pow((1-p),2));
       }


         if(*model==57) // poissongamma inflado 2 nuggets
       {
         //  Rprintf("%f %f %f %f  %f \n",nuis[0],nuis[1],nuis[2],nuis[3],nuis[4]);
           mui=exp(mean[i]);muj=exp(mean[j]);
           ai=mui*(1+mui/nuis[4]);aj=muj*(1+muj/nuis[4]);
           dd=sqrt(ai*aj)*corr_pois_gen((1-nuis[1])*corr,mui, muj, nuis[4]); // it's the pem  covariance
           p=pnorm(nuis[3],0,1,1,0);
           psj=pbnorm22(nuis[3],nuis[3],(1-nuis[2])*corr);
           p1=1-2*p+psj;
           rho[h]=p1*dd +  ai*aj*(p1-R_pow((1-p),2));
      
       }

           if(*model==58) // poissongamma inflado 1 nugget
       {
       // Rprintf("%f %f %f %f   \n",nuis[0],nuis[1],nuis[2],nuis[3]);
           mui=exp(mean[i]);muj=exp(mean[j]);
           ai=mui*(1+mui/nuis[2]);aj=muj*(1+muj/nuis[2]);
           dd=sqrt(ai*aj)*corr_pois_gen((1-nuis[0])*corr,mui, muj, nuis[2]); // it's the pem  covariance
           p=pnorm(nuis[1],0,1,1,0);
           psj=pbnorm22(nuis[1],nuis[1],(1-nuis[0])*corr);
           p1=1-2*p+psj;
           rho[h]=p1*dd +  ai*aj*(p1-R_pow((1-p),2));
      
       }
            h++;
        }}
    return;
}

// Computation of the correlations for spatial tapering:
void CorrelationMat_tap(double *rho,double *coordx, double *coordy,double *coordz, double *coordt, int *cormod,  double *nuis, double *par,double *radius,
  int *ns, int *NS)
{
  int i=0;// check the paramaters range:
  for(i=0;i<*npairs;i++) {
       rho[i]=CorFct(cormod,lags[i],0,par,0,0);}
  return;
}
// Computation of the correlations for kringing with  sparse matrix:
void CorrelationMat_dis_tap(double *rho,double *coordx, double *coordy, double *coordz, double *coordt, int *cormod,  double *nuis, double *par,double *radius,
  int *ns, int *NS, int *n1,int *n2, double *mu1,double *mu2,int  *model)
{
  int i=0;
  double p1,p2,psj,aux1,aux2,p,corr=0.0,dd,p00,p11;
  for(i=0;i<*npairs;i++) {

       corr=CorFct(cormod,lags[i],0,par,0,0);
  /***************************************************************/
  if(*model==2||*model==11||*model==14||*model==16||*model==45){
      psj=pbnorm22(mu1[i],mu2[i],(1-nuis[0])*corr);
      p1=pnorm(mu1[i],0,1,1,0);  p2=pnorm(mu2[i],0,1,1,0);

      if(*model==2||*model==11)       rho[i]=fmin2(n1[i],n2[i]) *(psj-p1*p2);
      if(*model==14)                  rho[i]=(psj-p1*p2)/((-psj+p1+p2)*p1*p2);
      if(*model==16)                  rho[i]=cov_binom_neg(n1[0],psj,p1,p2);
      if(*model==45)
         {
           p=pnorm(nuis[2],0,1,1,0);
           p00=pbnorm22(nuis[2],nuis[2],(1-nuis[1])*corr); p11=1-2*p+p00;
           dd=cov_binom_neg(n1[0],psj,p1,p2);
           rho[i]=p11*dd +  (n1[0]*n1[0]*(1-p1)*(1-p2)/(p1*p2)) * (p11-(1-p)*(1-p));
         }
    }
      /***************************************************************/
    if(*model==30||*model==36) // poisson
       {
          aux1=exp(mu1[i]);aux2=exp(mu2[i]);
          rho[i]=sqrt(aux1*aux2)*corr_pois((1-nuis[0])*corr,aux1, aux2);
       }

         if(*model==43||*model==44) // poisson inflado
       {
           aux1=exp(mu1[i]);aux2=exp(mu2[i]);
           dd=sqrt(aux1*aux2)*corr_pois((1-nuis[0])*corr,aux1, aux2); // it's the  covariance poisson

           p=pnorm(nuis[2],0,1,1,0);
           psj=pbnorm22(nuis[2],nuis[2],(1-nuis[1])*corr);
           p1=1-2*p+psj;
          rho[i]=p1*dd +   aux1*aux2*(p1-(1-p)*(1-p));
       }
        if(*model==46||*model==47) // Poisson gamma
       {
           aux1=exp(mu1[i]);aux2=exp(mu2[i]);
           p00=aux1*(1+aux1/nuis[1]);
           p11=aux2*(1+aux2/nuis[1]);
           rho[i]=sqrt(p00*p11)*corr_pois_gen((1-nuis[0])*corr,aux1, aux2, nuis[1]); // it's the  covariance
       }
   /***************************************************************/
      }
  return;
}

/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/****************** SPATIO-TEMPORAL CORRELATION MATRIX (upper trinagular) ***********************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/*************************************************************************************************/



// Computation of the correlations for spatio-temporal tapering:
void CorrelationMat_st_tap(double *rho,double *coordx, double *coordy, double *coordz, double *coordt, int *cormod, double *nuis,
  double *par,double *radius, int *ns, int *NS)
{
  int i=0;
  for(i=0;i<*npairs;i++) {
      rho[i]=CorFct(cormod,lags[i],lagt[i],par,0,0);
  }
return;
}

// Computation of the correlations for spatio-temporal tapering:
void CorrelationMat_st_dis_tap(double *rho,double *coordx, double *coordy, double *coordz, double *coordt, int *cormod,  double *nuis, double *par,double *radius,
  int *ns, int *NS, int *n1,int *n2, double *mu1,double *mu2,int  *model)
{
  int i=0;
  double dd,p,p1,p2,psj,aux1,aux2,p11,p00,corr=0.0;
  for(i=0;i<*npairs;i++) {

       corr=CorFct(cormod,lags[i],lagt[i],par,0,0);

  if(*model==2||*model==11||*model==14||*model==16||*model==45){
      p1=pnorm(mu1[i],0,1,1,0); p2=pnorm(mu2[i],0,1,1,0);
          psj=pbnorm22(mu1[i],mu2[i],(1-nuis[0])*corr);
      if(*model==2||*model==11)       rho[i]=fmin2(n1[i],n2[i])*(psj-p1*p1);
      if(*model==14)                  rho[i]=(psj-p1*p2)/((-psj+p1+p2)*p1*p2);
      if(*model==16)                  rho[i]=cov_binom_neg(n1[0],psj,p1,p2);
      if(*model==45)
         {
           p=pnorm(nuis[2],0,1,1,0);
           p00=pbnorm22(nuis[2],nuis[2],(1-nuis[1])*corr);p11=1-2*p+p00;
           dd=cov_binom_neg(n1[0],psj,p1,p2);
           rho[i]=p11*dd +  (n1[0]*n1[0]*(1-p1)*(1-p2)/(p1*p2)) * (p11-(1-p)*(1-p));
         }

}
    if(*model==30||*model==36) // poisson
       {

          aux1=exp(mu1[i]);   aux2=exp(mu2[i]);
          rho[i]=sqrt(aux1*aux2)*corr_pois((1-nuis[0])*corr,aux1, aux2);
       }
      if(*model==43||*model==44) // poisson inflated
       {
           aux1=exp(mu1[i]);aux2=exp(mu2[i]);
           dd=sqrt(aux1*aux2)*corr_pois((1-nuis[0])*corr,aux1, aux2); // it's the  covariance
           p=pnorm(nuis[2],0,1,1,0);
           psj=pbnorm22(nuis[2],nuis[2],(1-nuis[1])*corr);
           p1=1-2*p+psj;
           rho[i]=p1*dd +  aux1*aux2*(p1-(1-p)*(1-p));
       }

          if(*model==46||*model==47) // Poisson gamma
       {
           aux1=exp(mu1[i]);aux2=exp(mu2[i]);
           p00=aux1*(1+aux1/nuis[1]);
           p11=aux2*(1+aux2/nuis[1]);
           rho[i]=sqrt(p00*p11)*corr_pois_gen((1-nuis[0])*corr,aux1, aux2, nuis[1]); // it's the  covariance
       }


      }
  return;
}
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/****************** Dynamic SPATIO-TEMPORAL CORRELATION MATRIX (upper trinagular) ***************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/

// Computation of the upper (lower) triangular  correlation matrix: spatial-temporal case
void CorrelationMat_st_dyn2(double *rho, double *coordx, double *coordy,double *coordz,  double *coordt,int *cormod,
  double *nuis, double *par,double *radius, int *ns,int *NS)
{

  int i=0,j=0,t=0,v=0,h=0; double dd=0.0,tt=0.0;
  for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[v];j++){
           dd=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],
            coordz[(i+NS[t])],coordz[(j+NS[v])],*REARTH);
           rho[h]=CorFct(cormod,dd,0,par,t,v);
           h++;}}

    else {
         tt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          dd=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],
            coordz[(i+NS[t])],coordz[(j+NS[v])],*REARTH);

rho[h]=CorFct(cormod,dd,tt,par,t,v);
              h++;
              }}
    }}}
  return;
}

// Computation of the upper (lower) triangular  correlation matrix for discrete models: spatial-temporal case
void CorrelationMat_st_dyn_dis2(double *rho,double *coordx, double *coordy, double *coordz, double *coordt,  int *cormod,  double *mean,int *n,
  double *nuis, double *par,double *radius, int *ns, int *NS, int *model)

{
    int i=0,j=0,t=0,v=0,h=0; double dd=0.0,tt=0.0,corr=0.0;
       double psj,ai,aj,p1,p2,p,p11,p00,bi,bj;
for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[v];j++){
          dd=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],
            coordz[(i+NS[t])],coordz[(j+NS[v])],*REARTH);
               corr=CorFct(cormod,dd,0,par,0,0);

          if(*model==14||*model==16||*model==2||*model==11||*model==45){

                        ai=mean[i+ns[t]*t];aj=mean[j+ns[t]*v];
                        psj=pbnorm22(ai,aj,(1-nuis[0])*corr);
                        p1=pnorm(ai,0,1,1,0); p2=pnorm(aj,0,1,1,0);
      if(*model==2||*model==11)      
        { rho[h]=fmin2(n[i+NS[t]],n[(j+NS[v])])*(psj-p1*p2);          }               //binomial
      if(*model==14)                  rho[h]=(psj-p1*p2)/((-psj+p1+p2)*p1*p2);        //goemettric
      if(*model==16)                  rho[h]=cov_binom_neg(n[0],psj,p1,p2);         //binomialnegative
      if(*model==45)
      {
           p=pnorm(nuis[2],0,1,1,0);
           p00=pbnorm22(nuis[2],nuis[2],(1-nuis[1])*corr); p11=1-2*p+p00;
           dd=cov_binom_neg(n[0],psj,p1,p2);
           rho[h]=p11*dd +  (n[0]*n[0]*(1-p1)*(1-p2)/(p1*p2)) * (p11-(1-p)*(1-p));
      }
    }
if(*model==30||*model==36) {       //poisson
           ai=exp(mean[i+ns[t]*t]);
           aj=exp(mean[j+ns[t]*v]);
           rho[h]= sqrt(ai* aj)*corr_pois((1-nuis[0])*CorFct(cormod,dd,0,par,t,v),ai, aj);
         }

if(*model==46||*model==47) {       //poisson gamma
  
           bi=exp(mean[i+ns[t]*t]);
           bj=exp(mean[j+ns[t]*v]);
           ai=bi*(1+bi/nuis[1]);
           aj=bj*(1+bj/nuis[1]);
           rho[h]=sqrt(ai*aj)*corr_pois_gen((1-nuis[0])*CorFct(cormod,dd,0,par,t,v),bi, bj, nuis[1]);
         }

  if(*model==43||*model==44) //poisson inflado
       {
           ai=exp(mean[i+ns[t]*t]);
           aj=exp(mean[j+ns[t]*v]);
           dd=sqrt(ai*aj)*corr_pois((1-nuis[0])*corr,ai, aj); // it's the  covariance
           p=pnorm(nuis[2],0,1,1,0);
           psj=pbnorm22(nuis[2],nuis[2],(1-nuis[1])*corr);
           p1=1-2*p+psj;
          rho[h]=p1*dd+ai*aj*(p1-(1-p)*(1-p));
       }

  h++;}}
          else {
         tt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          dd=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],
            coordz[(i+NS[t])],coordz[(j+NS[v])],*REARTH);
          corr=CorFct(cormod,dd,tt,par,t,v);

            if(*model==14||*model==16||*model==2||*model==11||*model==45){
                         ai=mean[i+ns[v]*t];aj=mean[j+ns[v]*v];
                         psj=pbnorm22(ai,aj,(1-nuis[0])*corr);
                         p1=pnorm(ai,0,1,1,0); p2=pnorm(aj,0,1,1,0);

             if(*model==2||*model==11)       //rho[h]=n[0]*(psj-p1*p2);                    //binomial
             { rho[h]=fmin2(n[i+NS[t]],n[(j+NS[v])])*(psj-p1*p2);          } 
             if(*model==14)                  rho[h]=(psj-p1*p2)/((-psj+p1+p2)*p1*p2);    //goemettric
             if(*model==16)                  rho[h]=cov_binom_neg(n[0],psj,p1,p2); ;      // negative binomial
             if(*model==45)  {
                      p=pnorm(nuis[2],0,1,1,0);
                      p00=pbnorm22(nuis[2],nuis[2],(1-nuis[1])*corr);
                      p11=1-2*p+p00;
                      dd=cov_binom_neg(n[0],psj,p1,p2);
                      rho[h]=p11*dd +  (n[0]*n[0]*(1-p1)*(1-p2)/(p1*p2)) * (p11-(1-p)*(1-p));
         }
                }
  if(*model==30||*model==36) {          //poisson
        ai=exp(mean[i+ns[v]*t]);
        aj=exp(mean[j+ns[v]*v]);
          rho[h]= sqrt(ai* aj)*corr_pois((1-nuis[0])*corr,ai, aj);


  }
    if(*model==46||*model==47) {          //poissongamma

        double bi=exp(mean[i+ns[v]*t]);
        double bj=exp(mean[j+ns[v]*v]);
        ai=bi*(1+bi/nuis[2]);aj=bj*(1+bj/nuis[2]);
        rho[h]=sqrt(ai*aj)*corr_pois_gen((1-nuis[0])*corr,bi, bj, nuis[2]);

  }

      if(*model==43||*model==44)  //poisson inflado
       {
             ai=exp(mean[i+ns[v]*t]);
             aj=exp(mean[j+ns[v]*v]);
           dd=sqrt(ai*aj)*corr_pois((1-nuis[0])*corr,ai, aj); // it's the  covariance
           p=pnorm(nuis[2],0,1,1,0);
           psj=pbnorm22(nuis[2],nuis[2],(1-nuis[1])*corr);
           p1=1-2*p+psj;
          rho[h]=p1*dd+ai*aj*(p1-(1-p)*(1-p));
       }

 h++;}}
            }}}
    return;
}

/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/****************** SPATIO-DYNAMICAL BIVARIATE CORRELATION MATRIX (upper trinagular) **********/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/

// Computation of the upper (lower) triangular covariance matrix: bivariate case
void CorrelationMat_biv_dyn2(double *rho,double *coordx, double *coordy,double *coordz, double *coordt, int *cormod,  double *nuis,
  double *par,double *radius, int *ns,int *NS)
{
  int i=0,j=0,t=0,v=0,h=0;double dd=0.0;
    for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i;j<ns[t];j++){
          dd=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],
            coordz[(i+NS[t])],coordz[(j+NS[v])],*REARTH);
                rho[h]=CorFct(cormod,dd,0,par,t,v);
                h++;}}
    else {
         for(j=0;j<ns[v];j++){
          dd=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],
            coordz[(i+NS[t])],coordz[(j+NS[v])],*REARTH);
                rho[h]=CorFct(cormod,dd,0,par,t,v);
                h++;}}
    }}}
  return;
}
// Computation of the upper (lower) triangular covariance matrix: bivariate case
void CorrelationMat_biv_skew_dyn2(double *rho,double *coordx, double *coordy,double *coordz, double *coordt, int *cormod,
 double *nuis, double *par,double *radius, int *ns,int *NS)
{
 int i=0,j=0,t=0,v=0,h=0;
  double cc=0.0,lags=0.0;
  int N=2;

 // if(CheckCor(cormod,par)==-2){rho[0]=-2;return;}
  //  if(par[0]<=0 || par[1]<=0 || par[2]<0 || par[3]<0){rho[0]=-2;return;
     double *vari;
    vari=(double *) R_Calloc(N,double);
     double *sk;
    sk=(double *) R_Calloc(N,double);
    vari[0]=par[0]; vari[1]=par[1];
    par[0]=1;par[1]=1;
    sk[0]=nuis[2];sk[1]=nuis[3];
     //if(CheckCor(cormod,par)==-2){rho[0]=-2;return;}
    //if(par[0]<=0 || par[1]<=0 || par[2]<0 || par[3]<0){rho[0]=-2;return;}
           for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i;j<ns[t];j++){
            lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],
                coordz[(i+NS[t])],coordz[(j+NS[v])],*REARTH);
                cc=CorFct(cormod,lags,0,par,t,v);
   rho[h]=2*sk[t]*sk[v]*(sqrt(1-cc*cc)+cc*asin(cc)-1)/M_PI+sqrt(vari[t])*sqrt(vari[v])*cc;
            h++;}}
   else {
        for(j=0;j<ns[v];j++){
        lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],
            coordz[(i+NS[t])],coordz[(j+NS[v])],*REARTH);
                cc=CorFct(cormod,lags,0,par,t,v);
   rho[h]=2*sk[t]*sk[v]*(sqrt(1-cc*cc)+cc*asin(cc)-1)/M_PI+sqrt(vari[t])*sqrt(vari[v])*cc;
                h++;}}
    }}}
  return;
}
// Computation of the correlations for bivariate tapering:
void CorrelationMat_biv_tap(double *rho, double *coordx, double *coordy, double *coordz, double *coordt,int *cormod,
 double *nuis, double *par,double *radius, int *ns,int *NS)
{
  int i=0;

    for(i=0;i<*npairs;i++) {
      rho[i]=CorFct(cormod,lags[i],0,par,first[i],second[i]);}
  return;
}












/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/****************** SPATIO-Temporal(BIVARIATE) CORRELATION Vector *******************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/

//computation of correlation between a points and a vector (for kriging)
void Corr_c(double *cc,double *coordx, double *coordy,double *coordz, double *coordt, int *cormod, int *grid, double *locx,  double *locy,double *locz,
  int *ncoord, int *nloc,int *tloc,int *ns,int *NS,int *ntime, double *par, int *spt, int *biv, double *time,int *type, int *which,double *radius)
{

if(!spt[0]&&!biv[0])    {   //spatial case
int i=0,j=0,h=0;
double dis=0.0;
     for(j=0;j<(*nloc);j++){
         for(i=0;i<(*ncoord);i++){
           dis=dist(type[0],coordx[i],locx[j],coordy[i],locy[j],coordz[i],locz[j],radius[0]);
          cc[h]=CorFct(cormod,dis,0,par,0,0);
          h++;
          }}}

if(spt[0]) {
  int t=0,v=0,i=0,j=0,h=0;
double dis=0.0, dit=0.0;
        for(j=0;j<(*nloc);j++){
        for(v=0;v<(*tloc);v++){
           for(t=0;t<*ntime;t++){
                      dit=fabs(coordt[t]-time[v]);
                for(i=0;i<ns[t];i++){
                   //dis=dist(type[0],coordx[(i+NS[t])],locx[j],coordy[(i+NS[t])],locy[j],radius[0]);
                   dis=dist(type[0],coordx[(i+NS[t])],locx[j],coordy[(i+NS[t])],locy[j],coordz[(i+NS[t])],locz[j],radius[0]);
            
       cc[h]=CorFct(cormod,dis,dit,par,0,0);
        h++;}}}}
}
if(biv[0]) {
  int t,i,j,h=0;
double dis=0.0;
       for(j=0;j<(*nloc);j++){
            for(t=0;t<*ntime;t++){
                    for(i=0;i<ns[t];i++){
                     // dis=dist(type[0],coordx[(i+NS[t])],locx[j],coordy[(i+NS[t])],locy[j],radius[0]);
                                  dis=dist(type[0],coordx[(i+NS[t])],locx[j],coordy[(i+NS[t])],locy[j],coordz[(i+NS[t])],locz[j],radius[0]);
                                cc[h]=CorFct(cormod,dis,0,par,which[0],t);
                                h++;}}}}

}




///compute the covariance btwen loc to predict and locaton sites for binomial and geometric RF
void Corr_c_bin(double *cc,double *coordx, double *coordy,double *coordz, double *coordt, int *cormod, int *grid, double *locx,  double *locy, double *locz,int *ncoord, int *nloc,
                int *model,int *tloc,int *nn,int *n, int *ns,int *NS,int *ntime, double *mean,double *nuis, double *par, int *spt, int *biv, 
                double *time,int *type, int *which,double *radius, int *cop)
{


    if(!spt[0]&&!biv[0])  {   //spatial case
    int i=0,j=0,h=0;
    double ll,p,dd,dis=0.0,p1=0.0,p2=0.0,psj=0.0,ai=0.0,aj=0.0,corr=0.0,p00=0.0,p11=0.0,bi=0.0,bj=0.0;
            
for(j=0;j<(*nloc);j++){

                for(i=0;i<(*ncoord);i++){

                     dis=dist(type[0],coordx[i],locx[j],coordy[i],locy[j],coordz[i],locz[j],radius[0]);
                     corr=CorFct(cormod,dis,0,par,0,0);
                     ll=(1-nuis[0])*corr;

       /*****************************************************************/
      if(*model==2||*model==11||*model==19||*model==14||*model==16||*model==45||*model==30||*model==36)
                {
                        ai=mean[i];aj=mean[j];
                        psj=pbnorm22(ai,aj,ll);
                        p1=pnorm(ai,0,1,1,0); p2=pnorm(aj,0,1,1,0);

                       // compute the covariance!
                   if(*model==2||*model==11||*model==19) { if(*cop==0) cc[h]=fmin2(nn[j],n[i])*(psj-p1*p2);
                                                           if(*cop==1) cc[h]=ll;
                                                         }

                   if(*model==14)       {   if(*cop==0)   cc[h]=(psj-p1*p2)/((-psj+p1+p2)*p1*p2);
                                            if(*cop==1)  cc[h]=ll;//calculate_covariance((1-nuis[0])*corr, *model,n[i],n[j],p1,p2);}
                                        }
                   if(*model==16)          {   if(*cop==0)    cc[h]=cov_binom_neg(n[0],psj,p1,p2);
                                               if(*cop==1)    cc[h]=ll;//calculate_covariance((1-nuis[0])*corr, *model,n[i],n[j],p1,p2);}
                                           }
                   if(*model==45)   {
                             p=pnorm(nuis[2],0,1,1,0);
                             p00=pbnorm22(nuis[2],nuis[2],(1-nuis[1])*corr);p11=1-2*p+p00;
                             dd=cov_binom_neg(n[0],psj,p1,p2);
                             if(*cop==0) cc[h]=p11*dd +  (n[0]*n[0]*(1-p1)*(1-p2)/(p1*p2)) * (p11-(1-p)*(1-p));
                             if(*cop==1){ cc[h]=ll;}
                                    }
                   if(*model==30||*model==36)   //poisson
                               {
                         ai=exp(mean[i]);
                         aj=exp(mean[j]);
                         if(*cop==0) cc[h]=sqrt(ai*aj)*corr_pois(ll,ai, aj);
                         if(*cop==1) cc[h]=ll;//cc[h]=calculate_covariance((1-nuis[0])*corr, *model,n[i],n[j],ai,aj);}
                               }
                if(*model==46||*model==47)   //poisson gamma
                        {
                        bi=exp(mean[i]);bj=exp(mean[j]);
                        ai= bi*(1+bi/nuis[1]); aj= bj*(1+bj/nuis[1]);
                        if(*cop==0){cc[h]=sqrt(ai*aj)*corr_pois_gen(ll,bi, bj, nuis[1]); }
                        if(*cop==1){ cc[h]=ll;} 
                       }

                 if(*model==43||*model==44) // poisson zip 
       {
           ai=exp(mean[i]);aj=exp(mean[j]);
              if(*cop==0){
                 dd=sqrt(ai*aj)*corr_pois((1-nuis[0])*corr,ai, aj); // it's the  covariance
                 p=pnorm(nuis[2],0,1,1,0);
                 psj=pbnorm22(nuis[2],nuis[2],(1-nuis[1])*corr);
                 p1=1-2*p+psj;
                 cc[h]=p1*dd +  ai*aj*(p1-(1-p)*(1-p));
             }
        if(*cop==1){ cc[h]=ll;}} 
       }
               /*****************************************************************/
                    h++;
        }}
      }
if(*spt) { // spacetime
      int i=0,j=0,h=0,t=0,v=0;
    double ll,p,dd,dit=0.0,dis=0.0,p1=0.0,p2=0.0,psj=0.0,p11=0.0,p00=0.0,ai=0.0,aj=0.0,corr=0.0,mui,muj;
        for(j=0;j<(*nloc);j++){
        for(v=0;v<(*tloc);v++){
           for(t=0;t<*ntime;t++){
                      dit=fabs(coordt[t]-time[v]);
               for(i=0;i<ns[t];i++){
                 //  dis=dist(type[0],coordx[(i+NS[t])],locx[j],coordy[(i+NS[t])],locy[j],radius[0]);
                       dis=dist(type[0],coordx[(i+NS[t])],locx[j],coordy[(i+NS[t])],locy[j],coordz[(i+NS[t])],locz[j],radius[0]);
                    corr=CorFct(cormod,dis,dit,par,t,v);
                       ll=(1-nuis[0])*corr;
          /*****************************************************************/
                 if(*model==2||*model==11||*model==19||*model==14||*model==16||*model==45||*model==30||*model==36)
                {
                             ai=mean[j+*nloc * v];
                             aj=mean[i+*nloc * t];
                             psj=pbnorm22(ai,aj,ll);
                             //   psj=pbnorm(cormod,dis,dit,ai,aj,nuis[0],nuis[1],par,0);
                                p1=pnorm(ai,0,1,1,0); p2=pnorm(aj,0,1,1,0);

                   if(*model==2||*model==11||*model==19) {
                              if(*cop==0){  cc[h]=fmin2(nn[j],n[i+NS[t]])*(psj-p1*p2);}
                             if(*cop==1){   cc[h]=ll;} 

                                      }//binomial



                   if(*model==14)           { 
                                      if(*cop==0){     cc[h]=(psj-p1*p2)/((-psj+p1+p2)*p1*p2); }
                                       if(*cop==1){   cc[h]=ll;} 


                   } // geometric


                   if(*model==16)    {      
                               if(*cop==0){   cc[h]=cov_binom_neg(n[0],psj,p1,p2);}
                                if(*cop==1){   cc[h]=ll;} 
 
               }

                   if(*model==45)       {
                                p=pnorm(nuis[2],0,1,1,0);
                                p00=pbnorm22(nuis[2],nuis[2],ll);
                                p11=1-2*p+p00;
                                dd=cov_binom_neg(n[0],psj,p1,p2);
                                cc[i]=p11*dd +  (n[0]*n[0]*(1-p1)*(1-p2)/(p1*p2)) * (p11-(1-p)*(1-p));
         }
                }
                if(*model==30||*model==36) // poisson
                {
                             ai=exp(mean[j+*nloc * v]);
                             aj=exp(mean[i+*nloc * t]);
                              if(*cop==0) cc[h]=sqrt(ai*aj)*corr_pois(ll,ai, aj);
                              if(*cop==1)   cc[h]=ll;

                }
                  if(*model==46||*model==47) // poissongamma
                {
                             mui=exp(mean[j+*nloc * v]);
                             muj=exp(mean[i+*nloc * t]);
                              ai= mui*(1+mui/nuis[1]); aj= muj*(1+muj/nuis[1]);
                             if(*cop==0)   cc[h]=sqrt(ai*aj)*corr_pois_gen(ll,mui, muj, nuis[1]);
                              if(*cop==1)   cc[h]=ll;

                }
         if(*model==43||*model==44) // poisson inflated
       {
           ai=exp(mean[j+*nloc * v]);
           aj=exp(mean[i+*nloc * t]);
           dd=sqrt(ai*aj)*corr_pois(ll,ai, aj); // it's the  covariance
           p=pnorm(nuis[2],0,1,1,0);
           psj=pbnorm22(nuis[2],nuis[2],(1-nuis[1])*corr);
           p1=1-2*p+psj;
             if(*cop==0)   cc[h]=p1*dd +  ai*aj*(p1-(1-p)*(1-p));
               if(*cop==1)   cc[h]=ll;
       }
                 /*****************************************************************/
                    h++;}}}}
          }
      /* if(*biv) {

                // in case of an irregular grid of coordinates:
                for(j=0;j<(*nloc);j++){
                    for(i=0;i<*ncoord;i++){
                        dis=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                        for(t=0;t<*ntime;t++){
                            cc[h]=CorFct(cormod,dis,0,par,which[0],t);
                            h++;}}}
              }*/
}




 void Corr_c_tap(double *cc,double *cc_tap,double *coordx, double *coordy, double *coordz, double *coordt, int *cormod, int *cormodtap, int *grid, double *locx,  double *locy,double *locz,
                 double *mxd,double *mxt, int *ncoord, int *nloc, int *ns,int *NS,int*tloc,int *ntime, double *par, int *spt, int *biv, double *time,int *type,int *which,double *radius)
{
int i,j,h=0,*modtap;
double *partap,dis=0.0;

modtap=(int *) R_Calloc(1,int);*modtap=*cormodtap+1;
if(!spt[0]&&!biv[0])    {   //spatial case

partap=(double *) R_Calloc(1,double);
partap[0]=*mxd; // compact support
// in case of an irregular grid of coordinates:
    for(j=0;j<(*nloc);j++){
      for(i=0;i<(*ncoord);i++){
          // dis=dist(type[0],coordx[i],locx[j],coordy[i],locy[j],radius[0]);
               dis=dist(type[0],coordx[i],locx[j],coordy[i],locy[j],coordz[i],locz[j],radius[0]);
          cc[h]=CorFct(cormod,dis,0,par,0,0);
          cc_tap[h]=cc[h]*CorFct(modtap,dis,0,partap,0,0);
          h++;
          }}

    R_Free(partap);
}
else{    //spatio temporal  case or bivariate case
int t,v;double dit=0.0;
if(*spt) {
    partap=(double *) R_Calloc(4,double);
    partap[0]=mxd[0]; // spatial compact support
    partap[2]=mxd[1]; // delta1 param for non separable taper
    partap[3]=mxd[2]; // delta1 param for non separable taper
    partap[1]=mxt[0]; // temporal compact support

 // in case of an irregular grid of coordinates:
              for(j=0;j<(*nloc);j++){
        for(v=0;v<(*tloc);v++){
           for(t=0;t<*ntime;t++){
                      dit=fabs(coordt[t]-time[v]);
               for(i=0;i<ns[t];i++){
                //  dis=dist(type[0],coordx[(i+NS[t])],locx[j],coordy[(i+NS[t])],locy[j],radius[0]);
                   dis=dist(type[0],coordx[(i+NS[t])],locx[j],coordy[(i+NS[t])],locy[j],coordz[(i+NS[t])],locz[j],radius[0]);
            cc[h]=CorFct(cormod,dis,dit,par,t,v);
             cc_tap[h]=cc[h]*CorFct(modtap,dis,dit,partap,t,v);
            h++;}}}}

    R_Free(partap);}
    if(*biv)  {

        partap=(double *) R_Calloc(3,double);
        partap[0]=mxd[0]; // spatial compact support
        partap[1]=mxd[1]; // temporal compact support
        partap[2]=mxd[2]; // temporal compact support
        partap[3]=mxd[3]; // colocated taper

             for(j=0;j<(*nloc);j++){
            for(t=0;t<*ntime;t++){
                    for(i=0;i<ns[t];i++){ //for(i=0;i<*ncoord;i++){
                        //dis=dist(type[0],coordx[i],locx[j],coordy[i],locy[j],radius[0]);
                          dis=dist(type[0],coordx[i],locx[j],coordy[i],locy[j],coordz[i],locz[j],radius[0]);
                                cc[h]=CorFct(cormod,dis,0,par,which[0],t);
                                h++;}}}

        R_Free(partap);}
}
    R_Free(modtap);
}


// Derivatives with respect to R_power2 of the Cauchy correlation model:
double DCauchyPow(double R_power2, double scale, double rho)
{
  return -rho*log(R_pow(rho,-2/R_power2));
}
// Derivatives with respect to scale of the Cauchy correlation model:
double DCauchySc(double lag, double R_power2, double scale, double rho)
{
 return R_power2*rho*R_pow(rho, 2/R_power2)*R_pow(lag,2)/R_pow(scale,3);
}
// Derivatives with respect to scale of the Exponential correlation model:
double DExpoSc(double lag, double scale, double rho)
{
 return rho*lag/R_pow(scale,2);
}
// Derivatives with respect to scale of the Gaussian correlation model:
double DGaussSc(double lag, double scale, double rho)
{
  return 2*rho*R_pow(lag,2)/R_pow(scale,3);
}
// Derivatives with respect to R_power1 of the generalised Cauchy correlation model:
double DGenCauP1(double lag, double R_power1, double R_power2, double scale, double rho)
{
  if(lag)return R_power2*rho/R_power1*(log(1+R_pow(lag/scale,R_power1))/R_power1-
           R_pow(lag/scale,R_power1)*log(lag/scale)/
           (1+R_pow(lag/scale,R_power1)));
  else return 0.0;
}
// Derivatives with respect to R_power2 of the generalised Cauchy correlation model:
double DGenCauP2(double lag, double R_power1, double scale, double rho)
{
  return -rho*log(1+R_pow(lag/scale,R_power1))/R_power1;
}
// Derivatives with respect to scale of the generalised Cauchy correlation model:
double DGenCauSc(double lag, double R_power1, double R_power2, double scale, double rho)
{
  if(lag) return rho/(1+R_pow(lag/scale,2))*R_power2*R_pow(lag,R_power1)/R_pow(scale,R_power1+1);
  else return 0.0;
}
// Derivatives with respect to scale of the sferical correlation model:
double DSferiSc(double lag, double scale)
{
  if(lag<=scale) return 1.5*lag*(1-R_pow(lag/scale, 2))/R_pow(scale, 2);
  else return 0.0;
}

// Derivatives with respect to scale of the Wen1 correlation model:
double DWen1Sc(double lag, double scale, double smooth)
{

    if(lag<=scale) return (smooth+5)*(smooth+6)*lag*lag*R_pow(lag-scale,4)*R_pow(((scale-lag)/scale),smooth)/R_pow(scale,7);
    else return 0;
}

// Derivatives with respect to R_power of the Stable correlation model:
double DStabPow(double lag, double R_power, double scale, double rho)
{
  if(lag) return -rho*R_pow(lag/scale,R_power)*log(lag/scale);
  else return 0.0;
}
// Derivatives with respect to scale of the Stable correlation model:
double DStabSc(double lag, double R_power, double scale, double rho)
{
  if(lag) return rho*R_power*R_pow(lag/scale,R_power)/scale;
  else return 0.0;
}
// Derivatives with respect to scale of the Whittle-Matern correlation model:
double DWhMatSc(double eps, double lag, double scale, double smooth)
{
  if (lag){
    double pscale=0.0;
    pscale=(bessel_k(lag/(scale+eps),smooth,1)-
      bessel_k(lag/scale,smooth,1))/ eps;
    return R_pow(2,1-smooth)/gammafn(smooth)*R_pow(lag/scale,smooth)*
      (pscale-smooth*bessel_k(lag/scale,smooth,1)/scale);}
  else return 0;
}

// Derivatives with respect to scale of the wave model
double DWaveSc(double lag, double scale)
{
    if(lag==0) return 0.0;
    else return  sin(lag/scale)/lag-cos(lag/scale)/scale;
}

// Derivatives with respect to smooth of the Whittle-Matern correlation model:
double DWhMatSm(double eps, double lag, double scale, double smooth)
{
  if (lag){
    double psmooth=0.0;
    psmooth=(bessel_k(lag/scale,smooth+ eps,1)-
       bessel_k(lag/scale,smooth,1))/ eps;
    return R_pow(2,1-smooth)*R_pow(lag/scale,smooth)/gammafn(smooth)*
      ((log(lag/scale)-log(2)-digamma(smooth))*bessel_k(lag/scale,smooth,1)+psmooth);}
  else return 0;
}

// Derivatives with respect to smooth of the Wend1 correlation model:
double DWen1Sm(double lag, double scale, double smooth)
{
    if (lag<=scale) return R_pow(lag-scale,5)*R_pow((scale-lag)/scale,smooth)*( (log((scale-lag)/scale)*(smooth*lag+5*lag+scale))+lag)/R_pow(scale,6);
    else return 0;
}

double DMat_Cauchy_sc_t(double h,double u,double R_power2,double scale_s,double scale_t,double smooth)
{
  double arg=0,arg3=0;
  arg3=R_pow((1+R_pow(u/scale_t, 2)),-R_power2);
  if(h) arg=(R_pow(2, 1-smooth)/gammafn(smooth))*R_pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth, 1);
  else arg=1;
  return 2*R_pow(u,2)*R_power2*arg*arg3/(R_pow(scale_t,3)*(1+R_pow(u/scale_t, 2)));
}

double DMat_Cauchy_pw2(double h,double u,double R_power2,double scale_s,double scale_t,double smooth)
{
  double arg=0.0,arg2=0.0,arg3=0.0;
  arg3=R_pow((1+R_pow(u/scale_t, 2)),-R_power2);
  if(h) arg=(R_pow(2, 1-smooth)/gammafn(smooth))*R_pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth, 1);
  else  arg=1;
  if(1+R_pow(u/scale_t, 2)) arg2=arg*arg3*log(1+R_pow(u/scale_t, 2));
 return arg2;
}


double DMat_Cauchy_sc_s(double h,double u,double R_power2,double scale_s,double scale_t,double smooth)
{
  double arg1,arg2=0,arg3;
  arg3=R_pow((1+R_pow(u/scale_t, 2)),-R_power2);
  if(h) {  arg2=2*smooth*scale_s*bessel_k(h/scale_s, smooth,1)-h*bessel_k(h/scale_s, smooth+1,1);}
  arg1=R_pow(2, 1-smooth)*R_pow(h/scale_s, smooth)*arg3;
  return -arg1*arg2/(gammafn(smooth)*R_pow(scale_s,2));
}

double DMat_Cauchy_sm(double h,double u,double eps, double R_power2, double scale_s,double scale_t,double smooth)
{
  double arg=0.0,arg2=0.0,arg3=0.0,psmooth=0.0;
  arg3=R_pow((1+R_pow(u/scale_t, 2)),-R_power2);
  psmooth=(bessel_k(h/scale_s,smooth+ eps,1)-bessel_k(h/scale_s,smooth,1))/ eps;
  if(h) arg=(R_pow(2, 1-smooth)/gammafn(smooth))*R_pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth, 1);
  else arg=1;
  if(h) arg2=-arg3*arg*(log(2)+digamma(smooth)-log(h/scale_s)-psmooth/bessel_k(h/scale_s, smooth, 1));
  else arg2=0;
  return arg2;
}

double DMat_Exp_sc_t(double h,double u,double scale_s,double scale_t,double smooth)
{
  double arg=0;
  if(h) arg=(R_pow(2, 1-smooth)/gammafn(smooth))*R_pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth, 1);
 else arg=1;
  return arg*u*exp(-u/scale_t)/R_pow(scale_t,2);
}

double DMat_Exp_sc_s(double h,double u,double scale_s,double scale_t,double smooth)
{
double arg1=0,arg2=0;

 if(h){arg1=(R_pow(2, 1-smooth)/gammafn(smooth))*R_pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth, 1);
   arg2=h*(R_pow(2, 1-smooth)/gammafn(smooth))*R_pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth+1, 1);}
 else {arg1=1;
   arg2=0;}

 return -arg1*2*smooth*exp(-u/scale_t)/scale_s + arg2*exp(-u/scale_t)/R_pow(scale_s,2);
}

double DMat_Exp_sm(double h,double u,double eps,double scale_s,double scale_t,double smooth)
{
  double arg=0.0,psmooth=0.0,arg2=0.0;
  psmooth=(bessel_k(h/scale_s,smooth+ eps,1)-bessel_k(h/scale_s,smooth,1))/ eps;
  if(h){arg=(R_pow(2, 1-smooth)/gammafn(smooth))*R_pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth, 1);
    arg2=-exp(-u/scale_t)*arg*(log(2)+digamma(smooth)-log(h/scale_s)-psmooth/bessel_k(h/scale_s,smooth,1));}
  else arg2=0;
  return arg2;
}

/***************************************************/
/* derivative of bivariate wendland2 model */
/***************************************************/
double DWen1sep_biv_var1(double h,double var11,double var22,double nug11,double nug22, double scale, double col,double smoo,int c11,int c22)
{
    double rho=0.0;
    if((c11==0)&&(c22==0))                              rho= CorFunWend1_tap(h,scale,smoo);
    if((c11==0&&c22==1)||(c11==1&&c22==0))              rho=0.5*R_pow(var11,-0.5)*col*sqrt(var22)* CorFunWend1_tap(h,scale,smoo);
    return rho;
}
double DWen1sep_biv_var2(double h,double var11,double var22,double nug11,double nug22, double scale, double col,double smoo,int c11,int c22)
{
    double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0))             rho=0.5*R_pow(var22,-0.5)*col*sqrt(var11)* CorFunWend1_tap(h,scale,smoo);
    if((c11==1)&&(c22==1))                             rho= CorFunWend1_tap(h,scale,smoo);
    return rho;
}
double DWen1sep_biv_nug1(double h,double var11,double var22,double nug11,double nug22, double scale, double col,double smoo,int c11,int c22)
{
    double rho=0.0;;
    if((c11==0)&&(c22==0))  {if(h==0)   rho=1;}
    return rho;
}
double DWen1sep_biv_nug2(double h,double var11,double var22,double nug11,double nug22, double scale, double col,double smoo,int c11,int c22)
{
    double rho=0.0;;
    if((c11==1)&&(c22==1))  {if(h==0)   rho=1;}
    return rho;
}
double DWen1sep_biv_scale(double h,double var11,double var22,double nug11,double nug22, double scale, double col,double smoo,int c11,int c22)
{
    double rho=0.0;
    if((c11==0)&&(c22==0))                   rho=var11*DWen1Sc(h,scale,smoo);
    if((c11==0&&c22==1)||(c11==1&&c22==0))   rho=col*sqrt(var11*var22)*DWen1Sc(h,scale,smoo);
    if((c11==1)&&(c22==1))                   rho=var22*DWen1Sc(h, scale,smoo);
    return rho;
}
double DWen1sep_biv_col(double h,double var11,double var22,double nug11,double nug22, double scale, double col,double smoo,int c11,int c22)
{
    double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0))  rho=sqrt(var11*var22)* CorFunWend1_tap(h,scale,smoo);
    return rho;
}
double DWen1sep_biv_smoo(double h,double var11,double var22,double nug11,double nug22, double scale, double col,double smoo,int c11,int c22)
{
    double rho=0.0;
    if((c11==0)&&(c22==0))                   rho=var11*DWen1Sm(h,scale,smoo);
    if((c11==0&&c22==1)||(c11==1&&c22==0))   rho=col*sqrt(var11*var22)*DWen1Sm(h,scale,smoo);
    if((c11==1)&&(c22==1))                   rho=var22*DWen1Sm(h,scale,smoo);
    return rho;
}
/***************************************************/
/***************************************************/

/***************************************************/
/* derivative of bivariate matern separable  model */
/***************************************************/
double Dmatsep_biv_var1(double h,double var11,double var22,double nug11,double nug22, double scale,double smoo, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0)&&(c22==0))                  rho=CorFunWitMat(h, scale,smoo);
    if((c11==0&&c22==1)||(c11==1&&c22==0))  rho=0.5*R_pow(var11,-0.5)*col*sqrt(var22)*CorFunWitMat(h, scale,smoo);
    return rho;
}
double Dmatsep_biv_var2(double h,double var11,double var22,double nug11,double nug22, double scale,double smoo, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0))             rho=0.5*R_pow(var22,-0.5)*col*sqrt(var11)*CorFunWitMat(h, scale,smoo);;
    if((c11==1)&&(c22==1))                             rho=CorFunWitMat(h, scale,smoo);
    return rho;
}
double Dmatsep_biv_nug1(double h,double var11,double var22,double nug11,double nug22, double scale,double smoo, double col,int c11,int c22)
{
    double rho=0.0;;
    if((c11==0)&&(c22==0))  {if(h==0)   rho=1;}
    return rho;
}
double Dmatsep_biv_nug2(double h,double var11,double var22,double nug11,double nug22, double scale ,double smoo,double col,int c11,int c22)
{
    double rho=0.0;;
    if((c11==1)&&(c22==1))  {if(h==0)   rho=1;}
    return rho;
}
double Dmatsep_biv_scale(double h,double eps,double var11,double var22,double nug11,double nug22, double scale ,double smoo,double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0)&&(c22==0))                   rho=var11*DWhMatSc(eps,h,scale,smoo);
    if((c11==0&&c22==1)||(c11==1&&c22==0))   rho=col*sqrt(var11*var22)*DWhMatSc(eps,h,scale,smoo);
    if((c11==1)&&(c22==1))                   rho=var22*DWhMatSc(eps,h,scale,smoo);
    return rho;
}

double Dmatsep_biv_smo(double h,double eps,double var11,double var22,double nug11,double nug22, double scale,double smoo, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0)&&(c22==0))                   rho=var11*DWhMatSm(eps,h,scale,smoo);
    if((c11==0&&c22==1)||(c11==1&&c22==0))   rho=col*sqrt(var11*var22)*DWhMatSm(eps,h,scale,smoo);
    if((c11==1)&&(c22==1))                   rho=var22*DWhMatSm(eps,h,scale,smoo);
    return rho;
}
double Dmatsep_biv_col(double h,double var11,double var22,double nug11,double nug22, double scale,double smoo, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0))  rho=sqrt(var11*var22)*CorFunWitMat(h, scale,smoo);;
    return rho;
}
/***************************************************/
/***************************************************/
/***************************************************/
/* derivative of full bivariate matern model */
/***************************************************/
double DMat_biv_var1(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0)&&(c22==0))                     rho=CorFunWitMat(h, scale11,smoo11);
    if((c11==0&&c22==1)||(c11==1&&c22==0))     rho=0.5*R_pow(var11,-0.5)*col*sqrt(var22)*CorFunWitMat(h, scale12,smoo12);
    return rho;
}
double DMat_biv_var2(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0))             rho=0.5*R_pow(var22,-0.5)*col*sqrt(var11)*CorFunWitMat(h, scale12,smoo12);
    if((c11==1)&&(c22==1))                             rho=CorFunWitMat(h, scale22,smoo22);
    return rho;
}
double DMat_biv_nug1(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;;
    if((c11==0)&&(c22==0))  {if(h==0)   rho=1;}
    return rho;
}
double DMat_biv_nug2(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==1)&&(c22==1))  {if(h==0)   rho=1;}
    return rho;
}
double DMat_biv_scale1(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                       double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0)&&(c22==0))                 rho=var11*DWhMatSc(eps,h,scale11,smoo11);
    return rho;
}

double DMat_biv_scale2(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                       double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==1)&&(c22==1))                 rho=var22*DWhMatSc(eps,h,scale22,smoo22);
    return rho;
}
double DMat_biv_scale1_contr(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                             double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0)) rho=col*sqrt(var11*var22)*DWhMatSc(eps,h/0.5,scale11+scale22,0.5*(smoo11+smoo22));
    if((c11==0)&&(c22==0))                 rho=var11*DWhMatSc(eps,h,scale11,smoo11);
    return rho;
}
double DMat_biv_scale2_contr(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                             double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0)) rho=col*sqrt(var11*var22)*DWhMatSc(eps,h/0.5,scale11+scale22,0.5*(smoo11+smoo22));
    if((c11==1)&&(c22==1))                 rho=var22*DWhMatSc(eps,h,scale22,smoo22);
    return rho;
}
double DMat_biv_scale12(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                        double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0))    rho=col*sqrt(var11*var22)*DWhMatSc(eps,h,scale12,smoo12);
    return rho;
}
double DMat_biv_smo1(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0)&&(c22==0))                   rho=var11*DWhMatSm(eps,h,scale11,smoo11);
    return rho;
}
double DMat_biv_smo12(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                      double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0))   rho=col*sqrt(var11*var22)*DWhMatSm(eps,h,scale12,smoo12);
    return rho;
}
double DMat_biv_smo2(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==1)&&(c22==1))            rho=var22*DWhMatSm(eps,h,scale22,smoo22);
    return rho;
}
double DMat_biv_col(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                    double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0))  rho=sqrt(var11*var22)*CorFunWitMat(h, scale12,smoo12);
    return rho;
}


/***************************************************/
/* derivative of full bivariate wend (contr) model */
/***************************************************/
double DWen1_biv_var1(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                                             double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
  double rho=0.0;
    if((c11==0)&&(c22==0))                     rho=CorFunWend1_tap(h, scale11,smoo11);
    if((c11==0&&c22==1)||(c11==1&&c22==0))     rho=0.5*R_pow(var11,-0.5)*col*sqrt(var22)*CorFunWend1_tap(h, scale12,smoo12);
    return rho;
}
double DWen1_biv_var2(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
   double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0))             rho=0.5*R_pow(var22,-0.5)*col*sqrt(var11)*CorFunWend1_tap(h, scale12,smoo12);
    if((c11==1)&&(c22==1))                             rho=CorFunWend1_tap(h, scale22,smoo22);
    return rho;
}
double DWen1_biv_nug1(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;;
    if((c11==0)&&(c22==0))  {if(h==0)   rho=1;}
    return rho;
}
double DWen1_biv_nug2(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==1)&&(c22==1))  {if(h==0)   rho=1;}
    return rho;
}
double DWen1_biv_scale1(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                      double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0)&&(c22==0))                 rho=var11*DWen1Sc(h,scale11,smoo11);
    return rho;
}

double DWen1_biv_scale2(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                       double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
   double rho=0.0;
    if((c11==1)&&(c22==1))                 rho=var22*DWen1Sc(h,scale22,smoo22);
    return rho;

}

double DWen1_biv_scale1_contr(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                       double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0)&&(c22==0))                 rho=var11*DWen1Sc(h,scale11,smoo11);
    if((c11==0&&c22==1)||(c11==1&&c22==0)) rho=col*sqrt(var11*var22)*DWen1Sc(h,scale22+scale11,0.5*(smoo11+smoo12));
    return rho;
}


double DWen1_biv_scale2_contr(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                       double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0)) rho=col*sqrt(var11*var22)*DWen1Sc(h,scale22+scale11,0.5*(smoo11+smoo12));
    if((c11==1)&&(c22==1))                 rho=var22*DWen1Sc(h,scale22,smoo22);
    return rho;
}
double DWen1_biv_scale12(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                       double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0))    rho=col*sqrt(var11*var22)*DWen1Sc(h,scale12,smoo12);
    return rho;
}
double DWen1_biv_smo1(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                       double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0)&&(c22==0))                   rho=var11*DWen1Sm(h,scale11,smoo11);
    return rho;
}
double DWen1_biv_smo12(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0))   rho=col*sqrt(var11*var22)*DWen1Sm(h,scale12,smoo12);
    return rho;
}
double DWen1_biv_smo2(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
     double rho=0.0;
    if((c11==1)&&(c22==1))            rho=var22*DWen1Sm(h,scale22,smoo22);
    return rho;
}
double DWen1_biv_col(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                    double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0))  rho=sqrt(var11*var22)*CorFunWend1_tap(h,scale12,smoo12);
    return rho;
}

/***************************************************/
/* derivative of LMC  (contr) model */
/***************************************************/

double DLMC_contr_var1(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22,double col,int c11,int c22)
{
    double rho=0.0,smoo11;
    smoo11=CorFunStable(h, 1, scale11);
    if(h==0) {smoo11=smoo11+nug11;}
    if((c11==0)&&(c22==0))                  {rho=2*var11*smoo11;}
    if((c11==0&&c22==1)||(c11==1&&c22==0))  {rho=col*smoo11;}
    return rho;
}
double DLMC_contr_var2(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22,double col,int c11,int c22)
{
    double rho=0.0,smoo22;
    smoo22=CorFunStable(h, 1, scale22);
    if(h==0) {smoo22=smoo22+nug22;}
    if((c11==0&&c22==1)||(c11==1&&c22==0))  {rho=col*smoo22;}
    if((c11==1)&&(c22==1))                  {rho=2*var22*smoo22;}
    return rho;
}
double DLMC_contr_nug1(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22,double col,int c11,int c22)
{
    double rho=0.0,smoo11,smoo22;
    smoo11=CorFunStable(h, 1, scale11);
    smoo22=CorFunStable(h, 1, scale22);
    if(h==0) {smoo11=1;smoo22=0;}
    else   {smoo11=0;smoo22=0;}
    if((c11==0)&&(c22==0))                  {rho=R_pow(var11,2)*smoo11+R_pow(col,2)*smoo22;}
    if((c11==0&&c22==1)||(c11==1&&c22==0))  {rho=var11*col*smoo11+var22*col*smoo22;}
    if((c11==1)&&(c22==1))                  {rho=R_pow(col,2)*smoo11+R_pow(var22,2)*smoo22;}
    return rho;
}
double DLMC_contr_nug2(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22,double col,int c11,int c22)
{
    double rho=0.0, smoo11,smoo22;
    smoo11=CorFunStable(h, 1, scale11);
    smoo22=CorFunStable(h, 1, scale22);
    if(h==0) {smoo11=0;smoo22=1;}
     else   {smoo11=0;smoo22=0;}
    if((c11==0)&&(c22==0))                  {rho=R_pow(var11,2)*smoo11+R_pow(col,2)*smoo22;}
    if((c11==0&&c22==1)||(c11==1&&c22==0))  {rho=var11*col*smoo11+var22*col*smoo22;}
    if((c11==1)&&(c22==1))                  {rho=R_pow(col,2)*smoo11+R_pow(var22,2)*smoo22;}
    return rho;
}

double DLMC_contr_col(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22,double col,int c11,int c22)
{
    double rho=0.0, smoo11,smoo22;
    smoo11=CorFunStable(h, 1, scale11);
    smoo22=CorFunStable(h, 1, scale22);
    if(h==0) {smoo11=smoo11+nug11;smoo22=smoo22+nug22;}
    if((c11==0)&&(c22==0))                  {rho=2*col*smoo22;}
    if((c11==0&&c22==1)||(c11==1&&c22==0))  {rho=var11*smoo11+var22*smoo22;}
    if((c11==1)&&(c22==1))                  {rho=2*col*smoo11;}
    return rho;
}
double DLMC_contr_scale11(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22,double col,int c11,int c22)
{
    double rho=0.0, smoo11;
    smoo11=DStabSc(h, 1,scale11,CorFunStable(h, 1, scale11));
    if(h==0) {smoo11=smoo11+nug11;}
    if((c11==0)&&(c22==0))                  {rho=R_pow(var11,2)*smoo11;}
    if((c11==0&&c22==1)||(c11==1&&c22==0))  {rho=var11*col*smoo11;}
    if((c11==1)&&(c22==1))                  {rho=R_pow(col,2)*smoo11;}
    return rho;
}
double DLMC_contr_scale22(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22,double col,int c11,int c22)
{
    double rho=0.0,smoo22;
    smoo22=DStabSc(h, 1,scale22,CorFunStable(h, 1, scale22));
    if(h==0) {smoo22=smoo22+nug22;}
    if((c11==0)&&(c22==0))                  {rho=R_pow(col,2)*smoo22;}
    if((c11==0&&c22==1)||(c11==1&&c22==0))  {rho=var22*col*smoo22;}
    if((c11==1)&&(c22==1))                  {rho=R_pow(var22,2)*smoo22;}
    return rho;
}




void GradCorrFct(double rho, int *cormod, double eps, int *flag,
     double *grad, double h, double u, int c11, int c22,double *par)
{
  int i=0;
  double R_power=0.0, R_power1=0.0, R_power2=0.0, R_power_s=0, R_power_t=0;
  double scale=0.0, scale_s=0, scale_t=0, smooth=0.0;
  double var11=0.0, var22=0.0, nug11=0.0, nug22=0.0,scale11=0.0,scale12=0.0,scale22=0.0,col=0.0,smoo11=0.0,smoo12=0.0,smoo22=0.0;

  switch(*cormod)// Correlation functions are in alphabetical order
    {//spatial gradients of correlations:
    case 1:// Cauchy correlation function
      R_power2=par[0];
      scale=par[1];
      if(flag[0]==1){//R_power parameter
    grad[i]=DCauchyPow(R_power2,scale,rho);i++;}
      if(flag[1]==1)//scale parameter
  grad[i]=DCauchySc(h, R_power2, scale, rho);
      break;
    case 4:// Exponential correlation function
      scale=par[0];//scale parameter
      if(flag[0] == 1) grad[i]=DExpoSc(h, scale, rho);
      break;
    //case 6:// Gaussian correlation function
     // scale=par[0];//scale parameter
     // if(flag[0]==1) grad[i]=DGaussSc(h,scale,rho);
     // break;
    case 8:// Generalised Cuachy correlation function
      R_power1=par[0];
      R_power2=par[1];
      scale=par[2];
      if(flag[0]==1){//R_power1 parameter
  grad[i]=DGenCauP1(h,R_power1,R_power2,scale,rho);i++;}
      if(flag[1]==1){//R_power2 parameter
  grad[i]=DGenCauP2(h,R_power1,scale,rho);i++;}
      if(flag[2]==1)//scale parameter
  grad[i]=DGenCauSc(h,R_power1,R_power2,scale,rho);
      break;
    case 10:// Sferical correlation function
      scale=par[0];//scale parameter
      if(flag[0]==1) grad[i]=DSferiSc(h,scale);
      break;
   /* case 11:// wend2 correlation function
        scale=par[0];//scale parameter
        if(flag[0]==1) grad[i]=DwenSc(h,scale);
        break;*/
    case 12:// Stable correlation function
      R_power=par[0];
      scale=par[1];
      if(flag[0]==1){//R_power parameter
  grad[i]=DStabPow(h,R_power,scale,rho);i++;}
      //scale parameter
      if(flag[1]==1) grad[i]=DStabSc(h,R_power,scale,rho);
      break;
    case 14:// Whittle-Matern correlation function
      scale=par[0];
      smooth=par[1];
      if(flag[0]==1){//scale parameter
  grad[i]=DWhMatSc(eps,h,scale,smooth);i++;}
      //smooth parameter
      if(flag[1]==1) grad[i]=DWhMatSm(eps,h,scale,smooth);
      break;
   case 16:// wave  correlation function
      scale=par[0];//scale parameter
      if(flag[0]==1) grad[i]=DWaveSc(h,scale);
      break;
      //spatio-temproal gradients of correlations:

    case 84://Double Exponential
      scale_s=par[0];
      scale_t=par[1];
      if(flag[0]==1){//spatial-scale parameter
  grad[i]=DStabSc(h,1,scale_s,rho);i++;}
      //temporal-scale parameter
      if(flag[1]==1) grad[i]=DStabSc(u,1,scale_t,rho);
      break;
    /*case 86://Exponential-Gaussian
      scale_s=par[0];
      scale_t=par[1];
      smooth_s=par[2];
      smooth_t=par[3];
      // to do..../
      break;*/
    case 94://Stable-Stable
      R_power_s=par[0];
      R_power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      if(flag[0]==1){//spatial-R_power parameter
  grad[i]=DStabPow(h,R_power_s,scale_s,rho);i++;}
      if(flag[1]==1){//temporal-R_power parameter
  grad[i]=DStabPow(u,R_power_t,scale_t,rho);i++;}
      if(flag[2]==1){//spatial-scale parameter///
  grad[i]=DStabSc(h,R_power_s,scale_s,rho);i++;}
      if(flag[3]==1){//temporal-scale parameter
  grad[i]=DStabSc(u,R_power_t,scale_t,rho);}
     break;
    case 110:
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        scale_s=par[5];
        smooth=par[6];
        if(flag[0]==1){//first variance
            grad[i]=DWen1sep_biv_var1(h,var11,var22,nug11,nug22,scale_s,col,smooth,c11,c22);i++;}
        if(flag[1]==1){//second varuance
            grad[i]=DWen1sep_biv_var2(h,var11,var22,nug11,nug22,scale_s,col,smooth,c11,c22);i++;}
        if(flag[2]==1){//first nugget
            grad[i]=DWen1sep_biv_nug1(h,var11,var22,nug11,nug22,scale_s,col,smooth,c11,c22);i++;}
        if(flag[3]==1){//second nugget
            grad[i]=DWen1sep_biv_nug2(h,var11,var22,nug11,nug22,scale_s,col,smooth,c11,c22);i++;}
        if(flag[4]==1){//colocated
            grad[i]=DWen1sep_biv_col(h,var11,var22,nug11,nug22,scale_s,col,smooth,c11,c22);i++;}
        if(flag[5]==1){//scle
            grad[i]=DWen1sep_biv_scale(h,var11,var22,nug11,nug22,scale_s,col,smooth,c11,c22);i++;}
        if(flag[6]==1){//smooh
            grad[i]=DWen1sep_biv_smoo(h,var11,var22,nug11,nug22,scale_s,col,smooth,c11,c22);}
        break;
        case 114:    /// bi wen1 full
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        scale11=par[5];
        scale12=par[6];
        scale22=par[7];
        smoo11=par[8];
        smoo12=par[9];
        smoo22=par[10];
        if(flag[0]==1){//first variance
            grad[i]=DWen1_biv_var1(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[1]==1){//second varuance
            grad[i]=DWen1_biv_var2(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[2]==1){//first nugget
            grad[i]=DWen1_biv_nug1(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[3]==1){//second nugget
            grad[i]=DWen1_biv_nug2(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[4]==1){//colocated
            grad[i]=DWen1_biv_col(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[5]==1){//scale1
            grad[i]=DWen1_biv_scale1(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[6]==1){//scale12
            grad[i]=DWen1_biv_scale12(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[7]==1){//scale2
            grad[i]=DWen1_biv_scale2(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[8]==1){//smo1
            grad[i]=DWen1_biv_smo1(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[9]==1){//smo12
            grad[i]=DWen1_biv_smo12(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[10]==1){//smo2
            grad[i]=DWen1_biv_smo2(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);}
        break;
      case 120:  /// bi wen1 withcontr
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        scale11=par[5];
        scale22=par[6];
        smoo11=par[7];
        smoo22=par[8];
        if(flag[0]==1){//first variance
            grad[i]=DWen1_biv_var1(h,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[1]==1){//second varuance
            grad[i]=DWen1_biv_var2(h,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[2]==1){//first nugget
            grad[i]=DWen1_biv_nug1(h,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[3]==1){//second nugget
            grad[i]=DWen1_biv_nug2(h,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[4]==1){//colocated
            grad[i]=DWen1_biv_col(h,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[5]==1){//scale1
            grad[i]=DWen1_biv_scale1_contr(h,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[6]==1){//scale2
            grad[i]=DWen1_biv_scale2_contr(h,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[7]==1){//smo1
            grad[i]=DWen1_biv_smo1(h,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[8]==1){//smo2
            grad[i]=DWen1_biv_smo2(h,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);}
             break;
      case 122:  /// bi matern sep
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        scale=par[5];
        smooth=par[6];
        if(flag[0]==1){//first variance
            grad[i]=Dmatsep_biv_var1(h,var11,var22,nug11,nug22,scale,smooth,col,c11,c22);i++;}
        if(flag[1]==1){//second varuance
            grad[i]=Dmatsep_biv_var2(h,var11,var22,nug11,nug22,scale,smooth,col,c11,c22);i++;}
        if(flag[2]==1){//first nugget
            grad[i]=Dmatsep_biv_nug1(h,var11,var22,nug11,nug22,scale,smooth,col,c11,c22);i++;}
        if(flag[3]==1){//second nugget
            grad[i]=Dmatsep_biv_nug2(h,var11,var22,nug11,nug22,scale,smooth,col,c11,c22);i++;}
        if(flag[4]==1){//colocated
            grad[i]=Dmatsep_biv_col(h,var11,var22,nug11,nug22,scale,smooth,col,c11,c22);i++;}
        if(flag[5]==1){//scle
            grad[i]=Dmatsep_biv_scale(h,eps,var11,var22,nug11,nug22,scale,smooth,col,c11,c22);i++;}
        if(flag[6]==1){//smooth
            grad[i]=Dmatsep_biv_smo(h,eps,var11,var22,nug11,nug22,scale,smooth,col,c11,c22);}
        break;
      case 124:  /// LMC contr
        var11=par[0];
         col=par[1];
        var22=par[2];
        nug11=par[3];
        nug22=par[4];
        scale11=par[5];
        scale22=par[6];
        if(flag[0]==1){//first variance
            grad[i]=DLMC_contr_var1(h,eps,var11,var22,nug11,nug22,scale11,scale22,col,c11,c22);i++;}
        if(flag[2]==1){//second varuance
            grad[i]=DLMC_contr_var2(h,eps,var11,var22,nug11,nug22,scale11,scale22,col,c11,c22);i++;}
        if(flag[1]==1){//colocated
            grad[i]=DLMC_contr_col(h,eps,var11,var22,nug11,nug22,scale11,scale22,col,c11,c22);i++;}
        if(flag[3]==1){//first nugget
            grad[i]=DLMC_contr_nug1(h,eps,var11,var22,nug11,nug22,scale11,scale22,col,c11,c22);i++;}
        if(flag[4]==1){//second nugget
            grad[i]=DLMC_contr_nug2(h,eps,var11,var22,nug11,nug22,scale11,scale22,col,c11,c22);i++;}
        if(flag[5]==1){//scle1
            grad[i]=DLMC_contr_scale11(h,eps,var11,var22,nug11,nug22,scale11,scale22,col,c11,c22);i++;}
        if(flag[6]==1){//scale2
            grad[i]=DLMC_contr_scale22(h,eps,var11,var22,nug11,nug22,scale11,scale22,col,c11,c22);}
      break;
        case 118:   /// bi matern contr
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        scale11=par[5];
        scale22=par[6];
        smoo11=par[7];
        smoo22=par[8];
        if(flag[0]==1){//first variance
            grad[i]=DMat_biv_var1(h,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[1]==1){//second varuance
            grad[i]=DMat_biv_var2(h,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[2]==1){//first nugget
            grad[i]=DMat_biv_nug1(h,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[3]==1){//second nugget
            grad[i]=DMat_biv_nug2(h,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[4]==1){//colocated
            grad[i]=DMat_biv_col(h,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[5]==1){//scale1
            grad[i]=DMat_biv_scale1_contr(h,eps,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[6]==1){//scale2
            grad[i]=DMat_biv_scale2_contr(h,eps,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[7]==1){//smo1
            grad[i]=DMat_biv_smo1(h,eps,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
          if(flag[8]==1){//smo2
            grad[i]=DMat_biv_smo2(h,eps,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);}
        break;
        case 128:    /// bi matern full
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        scale11=par[5];
        scale12=par[6];
        scale22=par[7];
        smoo11=par[8];
        smoo12=par[9];
        smoo22=par[10];
        if(flag[0]==1){//first variance
            grad[i]=DMat_biv_var1(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[1]==1){//second varuance
            grad[i]=DMat_biv_var2(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[2]==1){//first nugget
            grad[i]=DMat_biv_nug1(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[3]==1){//second nugget
            grad[i]=DMat_biv_nug2(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[4]==1){//colocated
            grad[i]=DMat_biv_col(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[5]==1){//scale1
            grad[i]=DMat_biv_scale1(h,eps,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[6]==1){//scale12
            grad[i]=DMat_biv_scale12(h,eps,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[7]==1){//scale2
            grad[i]=DMat_biv_scale2(h,eps,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[8]==1){//smo1
            grad[i]=DMat_biv_smo1(h,eps,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[9]==1){//smo12
            grad[i]=DMat_biv_smo12(h,eps,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[10]==1){//smo2
            grad[i]=DMat_biv_smo2(h,eps,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);}
        break;
    }

}





// compute the gradient matrix (numcoord...) for the spatial field:
void DCorrelationMat_biv_tap(int *cormod,double *coordx, double *coordy, double *coordz, double *coordt,double *drho,double *eps,int *flagcor,int *nparcor,double *parcor,double *rho)
{
    int i=0,m=0,k=0;
    double *gradcor,*derho;
    //Initializes gradients:
    gradcor=(double *) R_alloc(*nparcor, sizeof(double));
     derho= (double *) R_Calloc(*npairs * *nparcor,double);
    k=0;
    for(i=0;i<*npairs;i++){
        GradCorrFct(rho[i],cormod,eps[0],flagcor,gradcor,lags[i],0,first[i],second[i],parcor);
        for(m=0;m<*nparcor;m++){derho[k]=gradcor[m];k++;}}
    k=0;
    for(i=0;i<*nparcor;i++){
        for(m=0;m<*npairs;m++){
            drho[k]=derho[i+m* *nparcor];k++;}}
     R_Free(derho);
    return;
}










void DCorrelationMat_biv(int *cormod,double *coordx, double *coordy,  double *coordz,double *coordt,double *drho,double *eps,int *flagcor,
      int *nparcor,double *parcor,double *rho)
{
 int i=0,j=0,h=0,k=0,s=0,t=0,v=0,st=0,npa=0;
 double *gradcor,*derho;

 st=ncoord[0] * *ntime;
 npa=0.5*st*(st-1)+st;
 gradcor=(double *) R_alloc(*nparcor, sizeof(double));
 derho= (double *) R_Calloc(npa * *nparcor,double);

for(i=0;i<ncoord[0];i++){
    for(t=0;t<*ntime;t++){
      for(j=i;j<ncoord[0];j++){
  if(i==j){
    for(v=t;v<*ntime;v++){
          GradCorrFct(rho[h],cormod,eps[0],flagcor,gradcor,0,0,t,v,parcor);
          h++;
   for(s=0;s<*nparcor;s++){
     derho[k]=gradcor[s];
     k++;}}}
     else {
          for(v=0;v<*ntime;v++){
               GradCorrFct(rho[h],cormod,eps[0],flagcor,gradcor,0,0,t,v,parcor);
               h++;
   for(s=0;s<*nparcor;s++){
     derho[k]=gradcor[s];
     k++;}}}}}}

 k=0;
 for(i=0;i<*nparcor;i++)
   for(j=0;j<npa;j++){
     drho[k]=derho[i+j* *nparcor];
     k++;}
    R_Free(derho);
 return;
}




void DCorrelationMat_biv2(int *cormod,double *coordx, double *coordy,  double *coordz,double *coordt,double *drho,double *eps,int *flagcor,
      int *nparcor,double *parcor,double *rho)
{
 int i=0,j=0,h=0,k=0,s=0,t=0,v=0,st=0,npa=0;
 double *gradcor,*derho,lags=0.0;

 st=ncoord[0] * *ntime;
 npa=0.5*st*(st-1)+st;
 gradcor=(double *) R_alloc(*nparcor, sizeof(double));
 derho= (double *) R_Calloc(npa * *nparcor,double);

for(i=0;i<ncoord[0];i++){
    for(t=0;t<*ntime;t++){
      for(j=i;j<ncoord[0];j++){
  if(i==j){
    for(v=t;v<*ntime;v++){
          GradCorrFct(rho[h],cormod,eps[0],flagcor,gradcor,0,0,t,v,parcor);
          h++;
   for(s=0;s<*nparcor;s++){
     derho[k]=gradcor[s];
     k++;}}}
     else {
        //  lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                    lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],coordz[i],coordz[j],*REARTH);
          for(v=0;v<*ntime;v++){
               GradCorrFct(rho[h],cormod,eps[0],flagcor,gradcor,lags,0,t,v,parcor);
               h++;
   for(s=0;s<*nparcor;s++){
     derho[k]=gradcor[s];
     k++;}}}}}}

 k=0;
 for(i=0;i<*nparcor;i++)
   for(j=0;j<npa;j++){
     drho[k]=derho[i+j* *nparcor];
     k++;}
    R_Free(derho);
 return;
}




// Computes the spatio-temporal variogram:
double Variogram(int *cormod, double h, double u, double nugget, double var, double *par)
{
  double vario=0.0;
  //Computes the variogram
  vario=var*(1-nugget)*(1-CorFct(cormod,h,u,par,0,0));
  return vario;
}


// Vector of spatio-temporal correlations:
void VectCorrelation(double *rho, int *cormod, double *h, int *nlags, int *nlagt,double *mean,int *model,
                     double *nuis,double *par, double *u,int *N)
{
  int i,j,t=0;
  double ai=0.0,aj=0.0,p1=0.0,p2=0.0,dd=0.0,psj=0.0,p11=0.0,p00=0.0,p=0.0,ccc=0.0;
  
  for(j=0;j<*nlagt;j++)
    for(i=0;i<*nlags;i++){
      if((*model==1)||(*model==10)||(*model==12)||(*model==21)||(*model==30)||(*model==36)||(*model==18)
        ||(*model==43)||(*model==44)||(*model==20)||(*model==33)||(*model==42)||(*model==46)||(*model==47)||
      (*model==22)||(*model==24)||(*model==26)||(*model==27)||(*model==29)||(*model==34)||(*model==38)
      ||(*model==39)||(*model==40)||(*model==41)||(*model==35)||(*model==37)||(*model==9)||(*model==50))


    {rho[t]=CorFct(cormod, h[i], u[j], par,0,0);}
           /***************************************/
      if(*model==2||*model==11||*model==19||*model==14||*model==16||*model==45)   // binomial  type I or  II or geometric case
      {
          ai=mean[i];aj=mean[j];
          ccc=CorFct(cormod, h[i], u[j], par,0,0);
          psj=pbnorm22(ai,aj,(1-nuis[0])*ccc);
          p1=pnorm(ai,0,1,1,0); p2=pnorm(aj,0,1,1,0);

          if(*model==2||*model==11||*model==19) {
            double denom = sqrt(p1*p2*(1-p1)*(1-p2));
            rho[t] = (denom != 0.0) ? (psj - p1*p2)/denom : 0.0;
        }
          if(*model==14)                        rho[t]= ((psj-p1*p2)/((-psj+p1+p2)*p1*p2))/(sqrt((1-p1)*(1-p2))/(p1*p2));  //covar12/sqrt(var1var2)
          if(*model==16)                        rho[t]=cov_binom_neg(N[0],psj,p1,p2)/(N[0]* sqrt((1-p1)*(1-p2))/(p1*p2));
          if(*model==45)
      {
           p=pnorm(nuis[2],0,1,1,0);
           p00=pbnorm22(nuis[2],nuis[2],(1-nuis[1])*ccc);
           p11=1-2*p+p00;
           dd=cov_binom_neg(N[0],psj,p1,p2);
           rho[t]=(p11*dd +  (N[0]*N[0]*(1-p1)*(1-p2)/(p1*p2)) * (p11-(1-p)*(1-p)))/
           (sqrt(N[0]*(1-p1)*(1-p)*(1+N[0]*p*(1-p1))*R_pow(p1,-2))*
            sqrt(N[0]*(1-p2)*(1-p)*(1+N[0]*p*(1-p2))*R_pow(p2,-2)));
      }
  }
      /***************************************/
      t++;}
  return;
}



// Vector of spatio-temporal correlations:


// Vector of spatio-temporal correlations:
void VectCorrelation_biv(double *rho, double *vario,int *cormod, double *h, int *nlags, int *nlagt,double *mean,int *model,
                     double *nuis,double *par, double *u,int *N)
{
    int i,j,p,t=0;
    for(j=0;j<2;j++)
    for(p=0;p<2;p++)
    for(i=0;i<*nlags;i++){
        rho[t]=CorFct(cormod, h[i],0, par,j,p);
        vario[t]=CorFct(cormod, 0,0, par,j,p)-CorFct(cormod, h[i],0, par,j,p);
        t++;}
    return;
}




//  Bounds for the bivariate matern
double bi_matern_bounds(double scale11,double scale22,double scale12,double nu11,double nu22,double nu12,double t,int c){
  double scale,l,inf;
  scale=R_pow(scale11,2*nu11)*R_pow(scale22,2*nu22)/  R_pow(scale12,4*nu12);
  if(!c) inf=R_pow(R_pow(scale12,2)+t*t,(2*nu12+2))/(R_pow((R_pow(scale11,2)+t*t),(nu11+1))*R_pow((R_pow(scale22,2)+t*t),(nu22+1)));
  else  inf=1;
  l=(gammafn(nu11+1)*gammafn(nu22+1)*R_pow(gammafn(nu12),2)*scale*inf)/(gammafn(nu11)*gammafn(nu22)*R_pow(gammafn(nu12+1),2));
  return(l);
}
