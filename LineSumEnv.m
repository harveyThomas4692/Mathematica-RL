(* ::Package:: *)

(* ::Title:: *)
(*LineSumEnv: RL Environment for Line Bundle Sums*)


(* ::Section::Closed:: *)
(*Start package*)


BeginPackage["LineSumEnv`"];


(* ::Section:: *)
(*Documentation*)


(* ::Subsection:: *)
(*General help*)


LineSumEnv::usage="\[FilledSmallSquare] LineSumEnv provides an environment for line bundle sums on Calabi-Yau three-folds in the context of heterotic model building.\n"<>
"\[FilledSmallSquare] The physics underlying this environment is explained in \!\(\*TemplateBox[{\"\\\"https://arxiv.org/pdf/1202.1757.pdf\\\"\",\"https://arxiv.org/pdf/1202.1757.pdf\"},\n\"HyperlinkURL\"]\) .\n"<>
"\[FilledSmallSquare] This environment plugs into RL or genetic algorithms.\n"<>
"\[FilledSmallSquare] A list of auxiliary modules is obtained by exectuing ?LineSumEnvAuxModules.\n"<>
"\[FilledSmallSquare] A list of main modules is obtained by executing ?LineSumEnvModules.\n"<>
"\[FilledSmallSquare] A list of global package options is obtained by executing ?LineSumEnvOptions.\n"<>
"\[FilledSmallSquare] The required data for a CY is loaded from a file. See ?LineSumEnvData and ?LoadCYData for more info.\n"<>
"\[FilledSmallSquare] States (= line bundle sums) are represented by integer matrices V, with each column the first Chern class of a line bundle.\n"<>
"\[FilledSmallSquare] Actions amount to adding +1 and -1 to two entries in the same row of a line bundle matrix V (so that c1(V)=0 is preserved).\n"<>
"\[FilledSmallSquare] Actions are represented by standard unit vectors, where the position of the 1 indicates the action.";

LineSumEnvAuxModules::usage="Auxiliary modules (execute ?modulename for more info):\n"<>
"\[FilledSmallSquare] MM[expr] converts matrices within a larger structure into matix form.\n"<>
"\[FilledSmallSquare] LoadCYData[name] loads the relevant data for a CY from a file.\n"<>
"\[FilledSmallSquare] RandomRow[rk,kmin,kmax], generator a random row of length rk which sums to zero.\n"<>
"\[FilledSmallSquare] RandomState[] generates a random state as a matrix.\n"<>
"\[FilledSmallSquare] Matrix2Bits[matrix,kmin,nbits] converts a matrix into a binary list.\n"<>
"\[FilledSmallSquare] Bits2Matrix[bitlist,kmin,nbits,ncol] converts a bitlist into a matrix.\n"<>
"\[FilledSmallSquare] RandomAction[] generates a random action, given by a standard unit vector.\n"<>
"\[FilledSmallSquare] BasicAct[statematrix,action] applies an action in standard unit vector form to a state in matrix form.\n"<>
"\[FilledSmallSquare] CheckSlope0[lbs,isec,kcone] checks if the slope 0 conditions for a line bundle sum have a solution in term of Kahler module t.\n"<>
"\[FilledSmallSquare] SCheckSlope0[lbs,skcone] checks if the slope 0 conditions for a line bundle sum have a solution in terms of dual Kahler moduli s.\n"<>
"\[FilledSmallSquare] TestSlope0[lbs,isec,range] checks a necessary condition for the slope 0 conditions to be satisfied.\n"<>
"\[FilledSmallSquare] ACohLine[linebundle,cicynum] computes line bundle cohomologies for select CICYs from formulae."; 

LineSumEnvModules::usage="Main modules (execute ?modulename for more info):\n"<>
"\[FilledSmallSquare] CompleteState[state,opt], computes all relevant quanties of a state list.\n"<>
"\[FilledSmallSquare] OrderState[state], brings a state into a canonical ordering.\n"<>
"\[FilledSmallSquare] LineSumInit[opt], generates a random state association.\n"<>
"\[FilledSmallSquare] LineSumAct[state,action,opt], acts an a state and computes new state and reward.";

LineSumEnvOptions::usage="Global package options:\n"<>
"\[FilledSmallSquare] \"SlopeZeroMeth\" -> max, if max=-1 (default) the full slope zero equations are solved in terms of the s variables (only possible if CY contains \"SKahlerCone\" data), if max=0 the full slope 0 equations are solved in terms of the t variables (only possible if CY containts \"KahlerCone\" data) and if max>0 the standard test for linear combintations with entries from -max,..,max is performed.\n"<>
"\[FilledSmallSquare] \"Rank\" -> integer, to specify the rank of the line bundle sum (default 5).\n"<>
"\[FilledSmallSquare] \"SymmOrder\" -> integer, to specify the order of the discrete symmetry (default 2).\n"<>
"\[FilledSmallSquare] \"MinEntry\" -> kmin, to specify the smallest integer allowed in a line bundle sums (default -3).\n"<>
"\[FilledSmallSquare] \"MaxEntry\" -> kmax, (default 4) which means the largest integer allowed in a line bundle sum is max(kmax,kmin-1+2^nbits).\n"<>
"\[FilledSmallSquare] \"NumBits\" -> nbits, to specify how many bits are used for one integer in a line bundle sum (default 3). The maximal entry in a line bundle sum is given by kmin-1+2^nbits.\n"<>
"\[FilledSmallSquare] \"ValueWeights\" -> list, to specify a list of 8 weights for the contributions to the intrinsic state value. The contributions are {anomaly,index,slope 0,equivariance,reduced structure group,non-zero h0 or h3,mirror generations,no Higgs pair,too many Higgs pairs}.\n"<>
"\[FilledSmallSquare] \"OutFormat\" -> string, to specify the output format for state associations. If set to \"Genetic\" or \"RL\" the items relevant for the respective algorithms are included, for \"Full\" (default) everything is included.\n"<>
"\[FilledSmallSquare] \"BoundaryPenalty\" -> real, to specify a (negative) number which is added to the reward if an action leads outside the range kmin,..,kmax.\n"<>
"\[FilledSmallSquare] \"TerminalBonus\" -> real, to specify a (positive) number added to the reward if an action leads to a terminal state.\n"<>
"\[FilledSmallSquare] \"RewardOffset\" -> real, to specify a (negative) number which provides the reward in case an action leads to a decrease in state value.\n"<>
"\[FilledSmallSquare] \"StepPenalty\" -> real, to specify a (negative) number which is added to every reward (in order to minimise episode lengths).\n"<>
"\[FilledSmallSquare] \"RewardPower\" -> p, to specify the power p in the computation of the reward so that reward=(newvalue-oldvalue)^p, provided newvalue>oldvalue.";

LineSumEnvData::usage="Data for CY manifolds is stored in files in $UserBaseDirectory<>\"Applications/LineSumEnvCYData/\". Data can be loaded in with the module LoadCYData and is stored in the variable CY, as an association with the following entries:\n"<>
"\[FilledSmallSquare] \"CicyNum\" -> integer, to specify the number for the underlying CY manifold.\n"<>
"\[FilledSmallSquare] \"NumPs\" -> integer, to specify the number of projective ambient spaces.\n"<>
"\[FilledSmallSquare] \"Conf\" -> matrix, to specify the CYs configuration matrix.\n"<>
"\[FilledSmallSquare] \"c2TX\" -> list, to specify the CYs second Chern class of the tangent bundle.\n"<>
"\[FilledSmallSquare] \"IntersecNumbers\" -> list, to specify the CY's triple intersection numbers.\n"<>
"\[FilledSmallSquare] \"KahlerCone\" -> list, to specify a list of vectors v, such that the Kahler cone is defined by v.t>0 for all v.\n"<>
"\[FilledSmallSquare] \"SKahlerCone\" -> list, so specify a list of vectors v, such that thet Kahler cone in the dual s variables is defined as v.s>0 for all v.\n"<>
"\[FilledSmallSquare] \"CYSymmetry\" -> list, to specify a list of permutation matrices which encode the CY symmetry (used for OrderState).";

CY::usage="Association which contains the relevant CY data. Load from file with LoadCYData. For entries see ?LineSumEnvData.";
LBScans::usage="Association which contains results of previous line bundle scans. Entries are of the form {cyis,symmorder}->list, where the list contains entries of the form {line bundle sum,number of mirror generations,number of Higgs pairs}.";


(* ::Subsection::Closed:: *)
(*Auxiliary Module help*)


MM::usage="MM[expr] converts matrices within a larger structure into matix form.";

LoadCYData::usage="LoadCYData[name] loads the relevant data for a CY from a file.\n"<>
"\[FilledSmallSquare] Loads data for CY with given name as an association into the variable CY.\n"<>
"\[FilledSmallSquare] The association keywords are listen under ?LineSumEnvData.\n"<>
"\[FilledSmallSquare] If the name \"LBScans\" is given then data from previous line bundle scans is loaded into variable LBScans (see ?LBScans for details).\n"<>
"\[FilledSmallSquare] The CY data file is loaded in from the directory $UserBaseDirectory<>\"Applications/LineSumEnvCYData/\".";

RandomRow::usage="RandomRow[rk,kmin,kmax], generator a random row of length rk which sums to zero, with entries in the range kmin to kmax and distributed as 1/(1+|k|).";

RandomState::usage="RandomState[] generates a random state as a matrix.";

Matrix2Bits::usage="Matrix2Bits[matrix,kmin,nbits] converts an integer matrix into a binary list, assuming minimal integer kmin and nbits bits per integer.";

Bits2Matrix::usage="Bits2Matrix[bitlist,kmin,nbits,ncol] converts a bitlist into an integer matrix, assuming minimal integer kmin, nbits bits per integer and a matrix with ncol columns.";

RandomAction::usage="RandomAction[] generates a random action as a standard unit vector.";

BasicAct::usage="BasicAct[state,action] applies an action in standard unit vector form to a state in matrix form.\n"<>
"\[FilledSmallSquare] Output is a list {newstate,isboundary}, where newstate is the new state in matrix form. If isboundary=False an action has been carried out, if isboundary=True the action would have crossed the boundary of the state space, so no action has been carried out and newstate equals the original state.";

CheckSlope0::usage="CheckSlope0[lbs,isec,kcone] checks if the slope 0 conditions for a line bundle sum have a solution in terms of Kahler moduli t.\n"<>
"\[FilledSmallSquare] lbs is the line bundle sum, isec the triple intersection numbers.\n"<>
"\[FilledSmallSquare] kcone defines the Kahler cone as the set of t with t.kcone[[a]]>0, for all a.\n"<>
"\[FilledSmallSquare] Output is zero if there is a solution or else equal to a number in (0,1] which measures how many maximial line bundle subsets have no common slope 0.";

SCheckSlope0::usage="SCheckSlope0[lbs,skcone] checks if the slope 0 conditions for a line bundle sum have a solution in term of dual Kahler moduli s.\n"<>
"\[FilledSmallSquare] lbs is the line bundle sum and skcone defines the Kahler cone in terms of the s variables as the set of s such that s.skcone[[a]]>0 for all a.\n"<>
"\[FilledSmallSquare] Output is zero if there is a solution or else equal to a number in (0,1] which measures how many maximial line bundle subsets have no common slope 0.";

TestSlope0::usage="TestSlope0[lbs,isec,range] checks a necessary condition for the slope 0 conditions to be satisfied.\n"<>
"\[FilledSmallSquare] lbs is the line bundle sum, isec the triple intersection numbers and linear combinations with entries from -range to range are considered.\n"<>
"\[FilledSmallSquare] Considers linear combintations of the matrices isec.Transpose[lbs][[a]] and checks if they have at least one positive and one negative entry.\n"<>
"\[FilledSmallSquare] This method is only applicable if the intersection numbers are >=0 and if the Kahler cone is the positive quadrant.\n"<>
"\[FilledSmallSquare] Output is a number in the range [0,1] which indicates which fraction of matrices fails the test.";

ACohLine::usage="ACohLine[linebundle,cicynum] computes line bundle cohomologies for select CICYs from formulae.\n"<>
"\[FilledSmallSquare] Here, linebundle is an integer vector representing the line bundle and cicynum is the number of the CICY. For linebundle={} a list of CICY numbers covered is returned.";


(* ::Subsection::Closed:: *)
(*Module help*)


CompleteState::usage="CompleteState[list,opt], computes all relevant quanties of a state.\n"<>
"\[FilledSmallSquare] The state can be given either as a state matix, a list of bits or a state association.\n"<>
"\[FilledSmallSquare] The output is a state association whose format is determined by the option \"OutFormat\" (see ?LineSumEnvOptions).\n"<>
"\[FilledSmallSquare] The global options \"SymmOder\", \"SlopeZeroMeth\" and \"OutFormat\" (see ?LineSumEnvOptions) can be overridded in the module call.";

OrderState::usage="OrderState[state], brings a state into a canonical ordering.";

LineSumInit::usage="LineSumInit[opt], generates a random state association.\n"<>
"\[FilledSmallSquare] Output is a state in association form, with format determined by the global option \"OutFormat\" (see ?LineSumEnvOptions) which can be overridded in the module call.";

LineSumAct::usage="LineSumAct[state,action,opt], acts on a state and computes new state and reward.\n"<>
"\[FilledSmallSquare] The state has to be provided in association form, the action as a standard unit vector.\n"<>
"\[FilledSmallSquare] The output is of the form <|\"State\"->state,\"NewState\"->newstate,\"Reward\"->reward,\"Terminal\"->terminal|>, where state and newstate are state associations.\n"<>
"\[FilledSmallSquare] The format of state and newstate is determined by the global option \"OutFormat\" (see ?LineSumEnvOptions) which can be overridded in the module call.";


(* ::Section:: *)
(*Initialize*)


(* ::Subsection::Closed:: *)
(*Start-up messages*)


Begin["`Private`"]
Print[Style["LineSumEnv: RL Environment for line bundle sums on CY manifolds",Underlined,FontColor -> Blue,TextAlignment -> Center, FontSize -> 14,FontFamily->"Times"]];
Print[Style["Execute \!\(\*
StyleBox[\"?\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"LineSumEnv\",\nFontWeight->\"Bold\"]\) for help.",FontColor -> Black,TextAlignment -> Center, FontSize -> 14,FontFamily->"Times"]];


(* ::Subsection::Closed:: *)
(*Options*)


Options[LineSumEnv]:={"SlopeZeroMeth"->-1,"Rank"->5,"SymmOrder"->2,"MinEntry"->-7,"MaxEntry"->8,"NumBits"->4,
                      "ValueWeights"->{1,1,1,1,1,1,1,1,1},"OutFormat"->"Full","BoundaryPenalty"->-1,
                      "TerminalBonus"->1,"RewardOffset"->-0.5,"StepPenalty"->0,"RewardPower"->1.2};
                      
LineSumEnv`CICYscovered={7890,7878,7889,7879,7861,7884,7887,7862,7844,7888,7880};

Print[Style["Default LineSumEnv option settings: ",FontColor -> Blue,TextAlignment -> Center, FontSize -> 14,FontFamily->"Times"],Options[LineSumEnv]];                                    


(* ::Section:: *)
(*Auxiliary Modules*)


(* ::Subsection::Closed:: *)
(*MM*)


MM[expr_]:=Module[{X},
 ReplaceAll[expr,X_/;MatrixQ[X]->MatrixForm[X,TableDepth->2]]]; 


(* ::Subsection::Closed:: *)
(*LoadCYData*)


LoadCYData[name_]:=Module[{dir,fn,flst},

dir=$UserBaseDirectory<>"/Applications/LineSumEnvCYData/";
fn=dir<>name;

If[name=="LBScans",
  LBScans=Get[fn];
  Print[Style["Line bundle data loaded into variable LBScans:",FontColor -> Blue,TextAlignment -> Center, FontSize -> 14,FontFamily->"Times"],Map[Length,LBScans]];
  Return[];
];

If[FileExistsQ[fn] && StringLength[name]>0,
  (* if file exist load in CY data *)
  CY=Get[fn];
  Print[Style["Data for "<>name<>" loaded into variable CY.",FontColor -> Blue,TextAlignment -> Left, FontSize -> 14,FontFamily->"Times"]],
  (* otherwise return list of available CYs *)
  flst=Map[Last,Map[FileNameSplit,FileNames[dir<>"*"]]];
  Print[Style["Available data: ",FontColor -> Blue,TextAlignment -> Center, FontSize -> 14,FontFamily->"Times"],flst];
]];

(* initialise CY data *)
LoadCYData["CICY7862"];           
LoadCYData[""];


(* ::Subsection::Closed:: *)
(*RandomRow*)


RandomRow[rk_,kmin_,kmax_]:=Module[{klst,kdis,row},

(* charge range and distribution *)
klst=Table[k,{k,kmin,kmax}];
kdis=Table[1/(1+Abs[k]^2),{k,kmin,kmax}];

row=RandomChoice[kdis->klst,rk];
While[Total[row]!=0,row=RandomChoice[kdis->klst,rk]];

row];


(* ::Subsection::Closed:: *)
(*RandomState*)


RandomState[]:=Module[{nP,rk,kmin,kmax},

(* options *)
nP=CY["NumPs"]; rk="Rank"/.Options[LineSumEnv];
kmin="MinEntry"/.Options[LineSumEnv]; 
kmax=Max["MaxEntry"/.Options[LineSumEnv],kmin+2^("NumBits"/.Options[LineSumEnv])-1];

Table[RandomRow[rk,kmin,kmax],{nP}]];


(* ::Subsection::Closed:: *)
(*Matrix2Bits*)


Matrix2Bits[matrix_,kmin_,nbits_]:=Module[{m},

m=Map[Most,matrix];
Flatten[Map[PadLeft[#,nbits]&,IntegerDigits[Flatten[m]-kmin,2]]]];


(* ::Subsection::Closed:: *)
(*Bits2Matrix*)


Bits2Matrix[bits_,kmin_,nbits_,ncol_]:=Module[{m},

m=Transpose[Partition[Map[FromDigits[#,2]&,Partition[bits,nbits]]+kmin,ncol-1]];
Transpose[Append[m,-Total[m]]]];


(* ::Subsection::Closed:: *)
(*RandomAction*)


RandomAction[]:=Module[{nP,rk,actdim},

nP=CY["NumPs"]; rk="Rank"/.Options[LineSumEnv];
actdim=nP*rk*(rk-1);

RandomChoice[IdentityMatrix[actdim]]];


(* ::Subsection::Closed:: *)
(*BasicAct*)


BasicAct[state_,action_]:=Module[
{nP,rk,kmin,kmax,pos,cols,row,colplus,colminus,newstate,entryplus,entryminus,boundary},

(* options *)
nP=CY["NumPs"]; rk="Rank"/.Options[LineSumEnv];
kmin="MinEntry"/.Options[LineSumEnv]; 
kmax=Max["MaxEntry"/.Options[LineSumEnv],kmin+2^("NumBits"/.Options[LineSumEnv])-1];

(* convert standard unit vector into row and two columns *)
pos=Position[action,1][[1,1]]-1;
cols=Quotient[pos,nP]; row=Mod[pos,nP]+1;
colminus=Quotient[cols,rk]+1; colplus=Mod[cols,rk]+1;
If[colminus>=colplus,colminus++];

(* compute new state *)
newstate=state; boundary=True;
entryplus=newstate[[row,colplus]]+1; entryminus=newstate[[row,colminus]]-1;

(* check if new state is outside the range *)
If[entryplus<=kmax && entryminus>=kmin,
  newstate[[row,colplus]]=entryplus;newstate[[row,colminus]]=entryminus; boundary=False];

{newstate,boundary}];


(* ::Subsection::Closed:: *)
(*CheckSlope0*)


CheckSlope0[lbs_,isec_,kcone_]:=Module[{tv,t,i,vT,nP,ineqs,vTred,eqs,slope0,slope0sub},

(* basic definitions *)
vT=Most[Transpose[lbs]]; nP=Length[First[vT]];
tv=Table[Subscript[t, i],{i,1,nP}];

(* inequalities which define Kahler cone *)
ineqs=Table[kcone[[a]] . tv>0,{a,1,Length[kcone]}];

(* row reduce line bundles for better performance *)
vTred=RowReduce[vT];

(* slope 0 equations *)
eqs=Table[isec . tv . tv . vTred[[a]]==0,{a,1,Length[vT]}];

(* check all equations for solutions *)
If[Length[FindInstance[Join[eqs,ineqs],tv]]>0,
  (* if there is a common solution, set slope measure to 0 *) 
  slope0=0,
  (* otherwise, check for subsets of equations *)
    slope0sub=Table[Length[FindInstance[Join[Delete[eqs,a],ineqs],tv]]>0,{a,1,Length[eqs]}];
    slope0=N[(1+Count[slope0sub,False])/(Length[vT]+1)]];

slope0];


(* ::Subsection::Closed:: *)
(*SCheckSlope0*)


SCheckSlope0[lbs_,skcone_]:=Module[{vT,nP,sv,s,i,ineqs,a,eqs,slope0,slope0sub},

(* basic definitions *)
vT=Most[Transpose[lbs]]; nP=Length[First[vT]];
sv=Table[Subscript[s, i],{i,1,nP}];

(* inequalities which define Kahler cone *)
ineqs=Table[skcone[[a]] . sv>0,{a,1,Length[skcone]}];

(* slope 0 equations *)
vT=RowReduce[vT];
eqs=Table[sv . vT[[a]]==0,{a,1,Length[vT]}];

(* check all equations for solutions *)
If[Reduce[Join[eqs,ineqs],sv]=!=False,
  (* if there is a common solution, set slope measure to 0 *) 
  slope0=0,
  (* otherwise, check for subsets of equations *)
    slope0sub=Table[Reduce[Join[Delete[eqs,a],ineqs],sv]=!=False,{a,1,Length[eqs]}];
    slope0=N[(1+Count[slope0sub,False])/(Length[vT]+1)]];

slope0];


(* ::Subsection::Closed:: *)
(*TestSlope0*)


TestSlope0[lbs_,isec_,range_]:=Module[{i,vT,coefflst,mlst,len,a,matlst,poslst,neglst,lst,slope0},

(* transpose and remove last line bundle, assuming c1=0 *)
vT=Most[Transpose[lbs]];len=Length[vT];

(* list of vectors for linear combinations *)
coefflst=Select[Tuples[Table[i,{i,-range,range}],len],Total[Abs[#]]>0&];

(* list of matrices whose linear combinations are considered *)
mlst=Table[isec . vT[[a]],{a,1,len}];

(* form linear combintations *)
matlst=Map[Flatten,Table[mlst . coefflst[[i]],{i,1,Length[coefflst]}]];

(* check if matrices have at least one positive and one negative entry *)
poslst=Apply[Or,Positive[matlst],{1}];neglst=Apply[Or,Negative[matlst],{1}];
lst=Thread[And[poslst,neglst]];

(* compute a measure for how well the test has gone *)    
slope0=N[Count[lst,False]/Length[lst]];

slope0];


(* ::Subsection::Closed:: *)
(*ACohLine *)


ACohLine[lb_,cicy_]:=Module[{},

If[Length[lb]==0,Return[LineSumEnv`CICYscovered]];

If[FreeQ[LineSumEnv`CICYscovered,cicy],Return[{}]];

If[cicy==7862,Return[CohLineTQ[lb]]];

Join[ACohLine01[lb,cicy],Reverse[ACohLine01[-lb,cicy]]]];

(* compute 0th and 1st cohomology from formulae *)
ACohLine01[lb_,cicynum_]:=Module[{lb1,h0,h1,k,k1,k2,k3,k4,ind,k0},

lb1=PadRight[lb,4];
k=lb1[[1]]; k1=k; k2=lb1[[2]]; k3=lb1[[3]]; k4=lb1[[4]];

Switch[cicynum,

7890,
     h0=Max[{KroneckerDelta[k,0]+5k^3/6+25k/6,0}];
     h1=0,

7878,
     h0=Max[{KroneckerDelta[k,0]+3k^3/2+9k/2,0}];
     h1=0,

7889,
     h0=Max[{KroneckerDelta[k,0]+4k^3/3+14k/3,0}];
     h1=0,

7879,
     h0=Max[{KroneckerDelta[k,0]+2k^3+5k,0}];
     h1=0,

7861,
     h0=Max[{KroneckerDelta[k,0]+8k^3/3+16k/3,0}];
     h1=0,

7884,
     ind=3(k1+k2)(2+k1 k2)/2;
     h0=Which[k1==0&&k2>=0,(1+k2)(2+k2)/2,
              k1>=0&&k2==0,(1+k1)(2+k1)/2,
              k1>0&&k2>0,ind,
              True,0];
     h1=Which[k1==0&&k2>0,(-1+k2)(-2+k2)/2,
              k1>0&&k2==0,(-1+k1)(-2+k1)/2,
              (k1<0 || k2<0)&&(k1+k2)>0,-ind,
               True,0],

7887,
     ind=(6k1(1+k2^2)+k2(11+k2^2))/3;
     h0=Which[k1>=0&&k2==0,k1+1,
              k1>=0&&k2>0,ind,
              k1<0&&k2==-4k1,-k1+1,
              k1<0&&k2>-4k1,32k1 (1-k1^2)/3+ind,
              True,0];
     h1=Which[k1<0&&k2==0,-k1-1,
              k1<-1&&k2>0&&-4k1>k2,-ind,
              k1<=-1&&k2==-4k1,-k1+1-ind,
              k1<=-1&&k2>-4k1,32k1 (1-k1^2)/3,
              True,0],

7844,
     ind=(6k1^2k2+9k1(1+k2^2)+k2(11+k2^2))/3;
     h0=Which[k1>=0&&k2==0,(1+k1)(2+k1)/2,
              k1>=0&&k2>0,ind,
              k1<0&&k2>=-6k1,8k1(2-3k1^2)+ind,
              True,0];
     h1=Which[k1>0&&k2==0,(1-k1)(2-k1)/2,
              k1>0&&(-k1+1)<k2&&k2<0,Max[{-ind,0}],
              k1<0&&(-k1-1)<=k2&&k2<-6k1&&k2>0,Max[{-ind,0}],
              k1<0&&k2>=-6k1,8k1(2-3k1^2),
              True,0],

7880,
     ind=(3k1 k2+k1 k3+k2 k3)k3+2(k1+k2)+3k3;
     If[k1>=k2,k0=k1;k1=k2;k2=k0];
     If[k3==0,
        h0=Which[k1>=0&&k2>=0,(1+k1)(1+k2),
                  True,0];
        h1=Which[k1>0&&k2>0,(-1+k1)(-1+k2),
                  k1<0&&k2>=0,(-1-k1)(1+k2),
                  True,0]];
     If[k3<0,
        h0=0;
        h1=Which[k1>0&&k2>0,Max[{-ind,0}],
                  True,0]];
     If[k3>0,
        h0=Which[k1==0&&k2==0,(1+k3)(2+k3)/2,
                  k1>=0&&k2>0,ind,
                  k1<0&&k2>=(-2k1-1)&&k3==-3k1,(1-k1)(1+2k1+k2),
                  k1<0&&k2>=(-2k1-1)&&k3>-3k1,9k1(1-k1^2)+ind,
                  True,0];
        h1=Which[k1==0&&k2==0,(-1+k3)(-2+k3)/2,
                  (k1<0&&k3<-3k1) || (k1<0&&k2<(-2k1-1)),Max[{-ind,0}],
                  k1<-1&&k2>=(-2k1-1)&&k3==-3k1,(1-k1)(1+2k1+k2)-ind,
                  k1<-1&&k2>=0&&k3!=-3k1,Max[{-ind,9k1(1-k1^2)}],
                  True,0]],

7888,
      ind=2 k1+(14 k2)/3+2 k1 k2^2+4/3 k2^3;
      h0=Which[k1>0&&k2>0,ind,
               k2>0&&k1+k2>0,ind+2k1/3-2k1^3/3,
               k1+k2<0||k2<0,0,
               k1>=0&&k2==0,k1+1,
               k2>0&&k1+k2==0,k2+1]; 
       h1="\!\(\*SuperscriptBox[\(h\), \(1\)]\)("<>ToString[lb]<>")",                         

_,
     h0="\!\(\*SuperscriptBox[\(h\), \(0\)]\)("<>ToString[lb]<>")";
     h1="\!\(\*SuperscriptBox[\(h\), \(1\)]\)("<>ToString[lb]<>")"];


{h0,h1}];

(* compute cohomology on tetra-quadric - probably needs work *)
CohLineTQ[lb_]:=Block[{h,n,p},
Do[Subscript[n, i]=Sort[lb][[i]],{i,1,4}];

h=Table[0,{4}];

Which[(Subscript[n, 1]<=-(4+2*(Subscript[n, 4]-2)))&&(Subscript[n, 2]==-(4+2*(Subscript[n, 4]-2)))&&(Subscript[n, 3]==-(4+2*(Subscript[n, 4]-2)))&&(Subscript[n, 4]>=2), 
p=Subscript[n, 4]-2;
h[[1]]=0;h[[2]]=0;h[[3]]=8*(p+1)*(p+2)*(p+3)+(p+1)*(-Subscript[n, 1]-(4+2*p)-1);
h[[4]]=Max[(Subscript[n, 4]+1),0]*\!\(
\*UnderoverscriptBox[\(\[Product]\), \(i = 1\), \(3\)]\(Max[\((\(-
\*SubscriptBox[\(n\), \(i\)]\) - 1)\)\ , 0]\)\)-(Max[(Subscript[n, 4]-1),0]*\!\(
\*UnderoverscriptBox[\(\[Product]\), \(i = 1\), \(3\)]\(Max[\((\(-
\*SubscriptBox[\(n\), \(i\)]\) + 1)\)\ , 0]\)\)-(8*(p+1)*(p+2)*(p+3)+(p+1)*(-Subscript[n, 1]-(4+2*p)-1))),
((Subscript[n, 1]<-(4+2*(Subscript[n, 4]-2)))&&(Subscript[n, 2]<-(4+2*(Subscript[n, 4]-2)))&&(Subscript[n, 3]<=-(4+2*(Subscript[n, 4]-2)))&&(Subscript[n, 4]>=2)\[Or](Subscript[n, 1]<=-(4+2*(Subscript[n, 4]-2))-2)&&(Subscript[n, 2]<=-(4+2*(Subscript[n, 4]-2))-2)&&(Subscript[n, 3]==-(4+2*(Subscript[n, 4]-2))+1)&&(Subscript[n, 4]>=2)),
p=Subscript[n, 4]-2;
h[[1]]=0;h[[2]]=0;h[[3]]=8*(p+1)*(p+2)*(p+3);
h[[4]]=Max[(Subscript[n, 4]+1),0]*\!\(
\*UnderoverscriptBox[\(\[Product]\), \(i = 1\), \(3\)]\(Max[\((\(-
\*SubscriptBox[\(n\), \(i\)]\) - 1)\)\ , 0]\)\)-(Max[(Subscript[n, 4]-1),0]*\!\(
\*UnderoverscriptBox[\(\[Product]\), \(i = 1\), \(3\)]\(Max[\((\(-
\*SubscriptBox[\(n\), \(i\)]\) + 1)\)\ , 0]\)\)-(8*(p+1)*(p+2)*(p+3))),
(Subscript[n, 1]<=-2)&&(Subscript[n, 2]==4+2*(-Subscript[n, 1]-2))&&(Subscript[n, 3]==4+2*(-Subscript[n, 1]-2))&&(Subscript[n, 4]>=4+2*(-Subscript[n, 1]-2)),
p=-Subscript[n, 1]-2;
h[[1]]=Max[(-Subscript[n, 1]+1),0]*\!\(
\*UnderoverscriptBox[\(\[Product]\), \(i = 2\), \(4\)]\(Max[\((
\*SubscriptBox[\(n\), \(i\)] - 1)\)\ , 0]\)\)-(Max[(-Subscript[n, 1]-1),0]*\!\(
\*UnderoverscriptBox[\(\[Product]\), \(i = 2\), \(4\)]\(Max[\((
\*SubscriptBox[\(n\), \(i\)] + 1)\)\ , 0]\)\)-(8*(p+1)*(p+2)*(p+3)+(p+1)*(Subscript[n, 4]-(4+2*p)-1)));
h[[2]]=8*(p+1)*(p+2)*(p+3)+(p+1)*(Subscript[n, 4]-(4+2*p)-1);h[[3]]=0;h[[4]]=0,
((Subscript[n, 1]<=-2)&&(Subscript[n, 2]>=4+2*(-Subscript[n, 1]-2))&&(Subscript[n, 3]>=4+2*(-Subscript[n, 1]-2))&&(Subscript[n, 4]>=4+2*(-Subscript[n, 1]-2))\[Or](Subscript[n, 1]<=-2)&&(Subscript[n, 2]==4+2*(-Subscript[n, 1]-2)-1)&&(Subscript[n, 3]>=4+2*(-Subscript[n, 1]-2)+2)&&(Subscript[n, 4]>=4+2*(-Subscript[n, 1]-2)+2)),
p=-Subscript[n, 1]-2;
h[[1]]=Max[(-Subscript[n, 1]+1),0]*\!\(
\*UnderoverscriptBox[\(\[Product]\), \(i = 2\), \(4\)]\(Max[\((
\*SubscriptBox[\(n\), \(i\)] - 1)\)\ , 0]\)\)-(Max[(-Subscript[n, 1]-1),0]*\!\(
\*UnderoverscriptBox[\(\[Product]\), \(i = 2\), \(4\)]\(Max[\((
\*SubscriptBox[\(n\), \(i\)] + 1)\)\ , 0]\)\)-(8*(p+1)*(p+2)*(p+3)));
h[[2]]=8*(p+1)*(p+2)*(p+3);h[[3]]=0;h[[4]]=0
]; 

If[h=={0,0,0,0},
h[[1]]=Max[\!\(
\*UnderoverscriptBox[\(\[Product]\), \(i = 1\), \(4\)]\(Max[\((
\*SubscriptBox[\(n\), \(i\)] + 1)\)\ , 0]\)\)-\!\(
\*UnderoverscriptBox[\(\[Product]\), \(i = 1\), \(4\)]\(Max[\((
\*SubscriptBox[\(n\), \(i\)] - 1)\)\ , 0]\)\),0]+Max[Max[(-Subscript[n, 1]+1),0]*\!\(
\*UnderoverscriptBox[\(\[Product]\), \(i = 2\), \(4\)]\(Max[\((
\*SubscriptBox[\(n\), \(i\)] - 1)\)\ , 0]\)\)-Max[(-Subscript[n, 1]-1),0]*\!\(
\*UnderoverscriptBox[\(\[Product]\), \(i = 2\), \(4\)]\(Max[\((
\*SubscriptBox[\(n\), \(i\)] + 1)\)\ , 0]\)\),0];
h[[2]]=Max[Max[(-Subscript[n, 1]-1),0]*\!\(
\*UnderoverscriptBox[\(\[Product]\), \(i = 2\), \(4\)]\(Max[\((
\*SubscriptBox[\(n\), \(i\)] + 1)\)\ , 0]\)\)- Max[(-Subscript[n, 1]+1),0]*\!\(
\*UnderoverscriptBox[\(\[Product]\), \(i = 2\), \(4\)]\(Max[\((
\*SubscriptBox[\(n\), \(i\)] - 1)\)\ , 0]\)\),0]+Max[Max[-Subscript[n, 1]+1,0]*Max[-Subscript[n, 2]+1,0]*Max[Subscript[n, 3]-1,0]*Max[Subscript[n, 4]-1,0]-Max[-Subscript[n, 1]-1,0]*Max[-Subscript[n, 2]-1,0]*Max[Subscript[n, 3]+1,0]*Max[Subscript[n, 4]+1,0],0];

h[[3]]=
Max[Max[-Subscript[n, 1]-1,0]*Max[-Subscript[n, 2]-1,0]*Max[Subscript[n, 3]+1,0]*Max[Subscript[n, 4]+1,0]-Max[-Subscript[n, 1]+1,0]*Max[-Subscript[n, 2]+1,0]*Max[Subscript[n, 3]-1,0]*Max[Subscript[n, 4]-1,0],0]+Max[Max[(Subscript[n, 4]-1),0]*\!\(
\*UnderoverscriptBox[\(\[Product]\), \(i = 1\), \(3\)]\(Max[\((\(-
\*SubscriptBox[\(n\), \(i\)]\) + 1)\)\ , 0]\)\)- Max[(Subscript[n, 4]+1),0]*\!\(
\*UnderoverscriptBox[\(\[Product]\), \(i = 1\), \(3\)]\(Max[\((\(-
\*SubscriptBox[\(n\), \(i\)]\) - 1)\)\ , 0]\)\),0];
h[[4]]=Max[Max[(Subscript[n, 4]+1),0]*\!\(
\*UnderoverscriptBox[\(\[Product]\), \(i = 1\), \(3\)]\(Max[\((\(-
\*SubscriptBox[\(n\), \(i\)]\) - 1)\)\ , 0]\)\)-Max[(Subscript[n, 4]-1),0]*\!\(
\*UnderoverscriptBox[\(\[Product]\), \(i = 1\), \(3\)]\(Max[\((\(-
\*SubscriptBox[\(n\), \(i\)]\) + 1)\)\ , 0]\)\),0]+Max[\!\(
\*UnderoverscriptBox[\(\[Product]\), \(i = 1\), \(4\)]\(Max[\((\(-
\*SubscriptBox[\(n\), \(i\)]\) + 1)\)\ , 0]\)\)-\!\(
\*UnderoverscriptBox[\(\[Product]\), \(i = 1\), \(4\)]\(Max[\((\(-
\*SubscriptBox[\(n\), \(i\)]\) - 1)\)\ , 0]\)\),0]];
h
];



(* ::Section:: *)
(*Modules*)


(* ::Subsection::Closed:: *)
(*CompleteState*)


Options[CompleteState]:={"OutFormat"->"Short"};

CompleteState[statelst_,opt___]:=Module[
{nP,rk,kmin,kmax,nbits,state,bits,stateass,isec,c2TX,vT,a,b,c2V,indLa,indV,slope0meth,kahlercone,i,tv,ineqs,eqs,
 slope0,slope0sub,coefflst,mlst,matlst,poslst,neglst,symmorder,vTg,equiv,cicynum,hLa,hV,hLaLb,h2V,hLaLbs,hVVs,lst,
 OX,nOX,statelst1,subc1,nsplits,valuelst,value,outformat,terminal,weights,t,nhiggs,vTred},

(* options *)
cicynum=CY["CicyNum"];nP=CY["NumPs"];isec=CY["IntersecNumbers"];
c2TX=CY["c2TX"]; kahlercone=CY["KahlerCone"];
rk="Rank"/.Options[LineSumEnv]; 
kmin="MinEntry"/.Options[LineSumEnv]; nbits="NumBits"/.Options[LineSumEnv];
kmax=Max["MaxEntry"/.Options[LineSumEnv],kmin-1+2^nbits];
symmorder="SymmOrder"/.{opt}/.Options[LineSumEnv];
slope0meth="SlopeZeroMeth"/.{opt}/.Options[LineSumEnv];
outformat="OutFormat"/.{opt}/.Options[LineSumEnv];
weights="ValueWeights"/.Options[LineSumEnv];

(* check if statelst is an association *)
If[AssociationQ[statelst],
  If[KeyExistsQ[statelst,"Input"],statelst1=statelst["Input"],statelst1=statelst["Bits"]],
  statelst1=statelst];

(* convert between matrix and bitlist form *)
If[MatrixQ[statelst1],state=statelst1; bits=Matrix2Bits[state,kmin,nbits]];
If[VectorQ[statelst],bits=statelst1; state=Bits2Matrix[bits,kmin,nbits,rk]];
vT=Transpose[state];
     
(* compute c2(V) *)
c2V=Sum[isec . vT[[a]] . vT[[a]],{a,1,rk}]/2;

(* compute index for each line bundle and total index *)
indLa=Table[a->(2*isec . vT[[a]] . vT[[a]] . vT[[a]]+c2TX . vT[[a]])/12,{a,1,rk}];
indV=Total[Map[Last,indLa]];

(* Check slope zero equations, as required by SlopeZeroMeth option *)
slope0=Switch[slope0meth,
-1,SCheckSlope0[state,CY["SKahlerCone"]],
 0,CheckSlope0[state,isec,kahlercone],
 _,TestSlope0[state,isec,slope0meth]];

(* check equivariance based on index *)
vTg=Gather[vT];
equiv=Sum[Mod[(2*isec . vTg[[a,1]] . vTg[[a,1]] . vTg[[a,1]]+c2TX . vTg[[a,1]])*Length[vTg[[a]]]/12,symmorder],{a,1,Length[vTg]}]/symmorder/Length[vTg];

(* count number of trivial line bundles in line bundle sum *)
OX=Table[0,{nP}]; nOX=Count[vT,OX];

(* check for non-trivial reduction of structure group *)
subc1=Map[Total,Subsets[vT,{2}]];
nsplits=Count[subc1,OX];

(* compute contributions to state value which do not involve cohomology *)
valuelst={
            10*Sum[Min[c2TX[[i]]-c2V[[i]],0],{i,1,nP}]/nP/kmax^2/rk,  (* penalty for violating the anomaly condition *)
            -400*Abs[indV+3*symmorder]/nP/kmax^3/rk,                (* penalty for wrong index *)
            -slope0/2,                                             (* penalty for violating slope 0 condition *)
            -equiv,                                                (* penalty for violating equivariance *)
            -(nOX+nsplits)/10};                                    (* penalty for reduced structure group *)

(* work out cohomologies if there is a cohomology code *)
If[MemberQ[LineSumEnv`CICYscovered,cicynum],
 hLa=Table[a->ACohLine[vT[[a]],cicynum],{a,1,rk}]; hV=Total[Map[Last,hLa]];
 hLaLb=Flatten[Table[{a,b}->ACohLine[vT[[a]]+vT[[b]],cicynum],{a,1,rk-1},{b,a+1,rk}],1];
 h2V=Total[Map[Last,hLaLb]]; nhiggs=h2V[[3]];
 hLaLbs=Flatten[Table[{a,b}->ACohLine[vT[[a]]-vT[[b]],cicynum],{a,1,rk},{b,1,rk}],1];
 hVVs=Total[Map[Last,hLaLbs]];
 
 (* compute contributions to state value which do involve cohomology*)
 valuelst=Join[valuelst,
              {
              -10*(hV[[1]]+hV[[4]]+h2V[[1]]+h2V[[4]])/nP/kmax^3/rk^2, (* penalty for non-zero h0 and h3 *)
              -10*hV[[3]]/nP/kmax^3,                            (* penalty for mirror generations *)
              (-1+HeavisideTheta[nhiggs-1/2])/10,              (* penalty for absence of Higgs pair *)
              -5*Max[nhiggs/symmorder-2,0]/nP/kmax^3/rk      (* penalty for too many Higgs pairs *)
              }],

  (* otherwise set cohomology entries to none *)
  hLa=None; hV=None; hLaLb=None; h2V=None;hLaLbs=None;hVVs=None;nhiggs=None;
];

(* multiply with weights and compute total value *)
valuelst=N[valuelst]*Take[weights,Length[valuelst]];
value=Total[valuelst]; 

(* determine if state is terminal *)
terminal=Abs[value]<10^(-5);

(* create output association for different formats *)          
stateass=Switch[outformat,

"Full",<|"Input"->state,"Bits"->bits,"c2(V)"->c2V,"Ind(La)"->indLa,"Ind(V)"->indV,"Slope0"->slope0,"Equiv"->equiv,
           "n(OX)"->nOX,"Splits"->nsplits,"h(La)"->hLa,"h(V)"->hV,"h(LaxLb)"->hLaLb,"h(\[CapitalLambda]2V)"->h2V,
           "h(LaxLb')"->hLaLbs,"h(VxV')"->hVVs,"ValueLst"->valuelst,"Value"->value,"Fitness"->value,"Terminal"->terminal|>,
           
"Genetic",<|"Input"->state,"Bits"->bits,"ValueLst"->valuelst,"Fitness"->value,"Terminal"->terminal|>,

"RL",<|"Input"->state,"ValueLst"->valuelst,"Value"->value,"Terminal"->terminal|>];
  
stateass];


(* ::Subsection::Closed:: *)
(*OrderState*)


OrderState[state_,opt___]:=Module[{cicysymm,st,kmin,nbits,rk,rowlst,lst,ordstate,outformat},

(* options *)
cicysymm=CY["CYSymmetry"];
kmin="MinEntry"/.Options[LineSumEnv]; nbits="NumBits"/.Options[LineSumEnv];
rk="Rank"/.Options[LineSumEnv];

(* prepare input *)
If[AssociationQ[state],st=state["Input"]];
If[MatrixQ[state],st=state];
If[VectorQ[state],st=Bits2Matrix[state,kmin,nbits,rk]];

(* generate all row permutations of state matrix encoded in cicysymm *)
rowlst=Map[Transpose,Map[Dot[#,st]&,cicysymm]];

(* generator all column permutation of these *)
lst=Map[Transpose,Flatten[Map[Permutations,rowlst],1]];

(* select state in standard form *)
ordstate=Sort[lst]//First;

CompleteState[ordstate,opt]];


(* ::Subsection::Closed:: *)
(*LineSumInit*)


LineSumInit[opt___]:=Module[{},

CompleteState[RandomState[],opt]];


(* ::Subsection::Closed:: *)
(*LineSumAct*)


LineSumAct[state_,action_,opt___]:=Module[{newstate,newstatemat,dvalue,reward,power,boundary,terminal},

(* options *)
power="RewardPower"/.Options[LineSumEnv];

(* compute new state matrix and complete to new state association *)
{newstatemat,boundary}=BasicAct[state["Input"],action];
newstate=CompleteState[newstatemat,opt];
terminal=newstate["Terminal"];

(* compute reward *)
dvalue=newstate["Value"]-state["Value"];
reward=Which[dvalue>0,dvalue^power,dvalue<=0,"RewardOffset"/.Options[LineSumEnv]]+("StepPenalty"/.Options[LineSumEnv]);

(* boundary penalty *)
If[boundary,reward=reward+("BoundaryPenalty"/.Options[LineSumEnv])];

(* terminal bonus *)
If[terminal,reward=reward+("TerminalBonus"/.Options[LineSumEnv])];

<|"State"->state,"NewState"->newstate,"Reward"->reward,"Terminal"->terminal|>];


(* ::Section::Closed:: *)
(*End Package*)


End[]
EndPackage[]
