(* ::Package:: *)

(* ::Title:: *)
(*LineBundleEnv: RL Environment for Line Bundles*)


(* ::Section::Closed:: *)
(*Start package*)


BeginPackage["LineBundleEnv`"];


(* ::Section::Closed:: *)
(*Documentation*)


(* ::Subsection::Closed:: *)
(*General help*)


LineBundleEnv::usage="\[FilledSmallSquare]  This package realises an RL environment for line bundles on Calabi-Yau manifolds.\n"<>
"\[FilledSmallSquare]  Basic parameters of the environment are defined in Options{LineBundleEnv], including the definition of the enrironment.\n"<>
"\[FilledSmallSquare]  Basic modules are listed under LineBundleEnvAuxModules and the main package modules under LineBundleEnvModules.\n"<>
"\[FilledSmallSquare]  The central module is the module LineBundleAct, which, given a state and an action, determines a new state and computes the reward.\n"<>
"\[FilledSmallSquare]  For the purpose of linking to an RL system a state is given in association form <|\"Input\"->statevector,\"Terminal\"->True/False,\"Value\"->value|>.\n"<>
"\[FilledSmallSquare]  An action (by LineBundleAct) leads to an output of the form <|\"State\"->state association,\"NewState\"->state association,\"Action\"->action vector,\"Reward\"->reward|>.";

LineBundleEnvAuxModules::usage="Auxiliary modules: \n"<>
"\[FilledSmallSquare]  RandomState[] generates a random state as a matrix of size nPx(rk(B)+rk(C)).\n"<>
"\[FilledSmallSquare]  RandomAction[] generates a random action as a unit vector.\n"<>
"\[FilledSmallSquare]  BasicAct[state,action] applies an action in vector form to a state in matrix form.\n"<>
"\[FilledSmallSquare]  Vector2Bits[vector,kmin,nbits] converts an integer vector into a list of bits.\n"<>
"\[FilledSmallSquare]  Bits2Vector[bits,kmin,nbits] converts a bit list into an integer vector.";

LineBundleEnvModules::usage="Modules: \n"<>
"\[FilledSmallSquare]  CompleteState[state] computes all properties of a state.\n"<>
(* "\[FilledSmallSquare]  OrderState[state], brings the line bundle into a canonical order.\n"<> *)
"\[FilledSmallSquare]  LineBundleAct[currentstate,action,opt] performs an action and determines the reward.\n"<>
"\[FilledSmallSquare]  LineBundleInit[] to produce an initial (random) state.";


(* ::Subsection::Closed:: *)
(*Auxiliary Module help*)


MM::usage="MM[list] converts matrices within the list into matix form.";

RandomState::usage="RandomState[] generates a random state as a vector.";

RandomAction::usage="RandomAction[] generates a random action as a unit vector.";

BasicAct::usage="BasicAct[state,action] applies an action in vector form to a state in vector form.";

Vector2Bits::usage="Vector2Bits[vector,kmin,nbits] converts an integer vector into a list of bits.";

Bits2Vector::usage="Bits2Vector[bits,kmin,nbits] converts a bit list into an integer vector.";


(* ::Subsection::Closed:: *)
(*Module help*)


CompleteState::usage="CompleteState[state] computes all properties of a line bundle state.\n"<>
"The state is given as a line bundle vector or an assocaciation <|\"Input\"->line bundle vector,...|>.\n"<>
"Output is a self-explanatory association with all relevant properties of the line bundle.";

LineBundleInit::usage="LineBundleInit[] to produce an initial (random) line bundle state.";

LineBundleAct::usage="LineBundleAct[currentstate,action,opt] performs an action and determines the reward.\n"<>
"\[FilledSmallSquare]  States are provided in association format <|\"Input\"->statevector,\"Terminal\"->True/False,\"Value\"->value|> and actions as a unit vector.\n"<>
"\[FilledSmallSquare]  The output is presented in the form <|\"State\"->stateassociation,\"NewState\"->stateassociation,\"Action\"->actionvector,\"Reward\"->reward\>.";

(* OrderState::usage="OrderState[state], brings the line bundles into a canonical order."; *)


(* ::Section::Closed:: *)
(*Initialize*)


(* ::Subsection::Closed:: *)
(*Start-up messages*)


Begin["`Private`"]
Print[Style["LineBundleEnv: RL Environment for line bundles on CY manifolds",Underlined,FontColor -> Blue,TextAlignment -> Center, FontSize -> 14,FontFamily->"Times"]];
Print[Style["Execute \!\(\*
StyleBox[\"?\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"LineBundleEnv\",\nFontWeight->\"Bold\"]\) for help. Current LineBundleEnv options:",FontColor -> Black,TextAlignment -> Center, FontSize -> 14,FontFamily->"Times"]];


(* ::Subsection::Closed:: *)
(*Options*)


Options[LineBundleEnv]:={"CicyNum"->7884,"NumPs"->2,"Conf"->{{3},{3}},"c2TX"->{36,36},"IntersecNumbers"->{{{0,3},{3,3}},{{3,3},{3,0}}},
                    "TargetIndex"->18,"EntryRange"->{-7,8},"BoundaryPenalty"->-1,"TerminalBonus"->2,"RewardOffset"->-1,"StepPenalty"->0,
                    "NumBits"->4};
Print[Style[Options[LineBundleEnv],"Output"]];              


(* ::Section::Closed:: *)
(*Auxiliary Modules*)


(* ::Subsection::Closed:: *)
(*RandomState*)


RandomState[]:=Module[{nP,range,klst,kdis,lb},

(* options *)
nP="NumPs"/.Options[LineBundleEnv]; range="EntryRange"/.Options[LineBundleEnv];

(* charge range and distribution *)
klst=Table[k,{k,range[[1]],range[[2]]}];
kdis=Table[1/(1+Abs[k]),{k,range[[1]],range[[2]]}];

(* generate random line bundle *)
lb=RandomChoice[kdis->klst,nP];

lb];


(* ::Subsection::Closed:: *)
(*RandomAction*)


RandomAction[]:=Module[{nP},

(* options *)
nP="NumPs"/.Options[LineBundleEnv];

RandomChoice[IdentityMatrix[2*nP]]];


(* ::Subsection::Closed:: *)
(*BasicAct*)


BasicAct[linebundle_,action_]:=Module[{actionvec,newlb,range,numbits,kmin,kmax},

(* options *)
range="EntryRange"/.Options[LineBundleEnv]; numbits="NumBits"/.Options[LineBundleEnv];
kmin=range[[1]]; kmax=Min[range[[2]],range[[1]]+2^numbits-1];

(* compute new line bundle *)
actionvec=Partition[action,2]/.{{0,0}->0,{1,0}->1,{0,1}->-1};
newlb=linebundle+actionvec;

(* check if new line bundle is outside bounaries *)
If[Max[newlb]>kmax || Min[newlb]<kmin,Return[linebundle],Return[newlb]]];


(* ::Subsection::Closed:: *)
(*Vector2Bits*)


Vector2Bits[vector_,kmin_,nbits_]:=Flatten[Map[PadLeft[#,nbits]&,IntegerDigits[vector-kmin,2]]];


(* ::Subsection::Closed:: *)
(*Bits2Vector*)


Bits2Vector[bits_,kmin_,nbits_]:=Map[FromDigits[#,2]&,Partition[bits,nbits]]+kmin;


(* ::Section::Closed:: *)
(*Modules*)


(* ::Subsection::Closed:: *)
(*CompleteState*)


CompleteState[state_,opt___]:=Module[{nP,isec,c2TX,targetind,lb,ind,range,maxentry,value,numbits,bits},

(* options *)
nP="NumPs"/.Options[LineBundleEnv]; isec="IntersecNumbers"/.Options[LineBundleEnv];
c2TX="c2TX"/.Options[LineBundleEnv]; targetind="TargetIndex"/.Options[LineBundleEnv];
range="EntryRange"/.Options[LineBundleEnv]; maxentry=Max[Abs[range]];
numbits="NumBits"/.Options[LineBundleEnv];

(* extract line bundle *)
If[AssociationQ[state],lb=state["Input"],lb=state];
If[Length[state]==nP,
  lb=state; bits=Vector2Bits[state,range[[1]],numbits],
  bits=state; lb=Bits2Vector[state,range[[1]],numbits]];

(* compute line bundle index *)
ind=isec.lb.lb.lb/6+c2TX.lb/12;

(* compute value *)
value=-10*Abs[ind-targetind]/maxentry^3/nP^3//N;

<|"Input"->lb,"Bits"->bits,"ind"->ind,"Value"->value,"Fitness"->value,"Terminal"->(value==0)|>];


(* ::Subsection::Closed:: *)
(*LineBundleInit*)


LineBundleInit[opt___]:=Module[{lb,cstate},

lb=RandomState[];
cstate=CompleteState[lb];

cstate];


(* ::Subsection::Closed:: *)
(*LineBundleAct*)


LineBundleAct[state_,action_,opt___]:=Module[{range,lb,newlb,newstatecompl,newstate,dvalue,reward},

(* options *)
range="EntryRange"/.Options[LineBundleEnv];

(* determine new state *)
lb=state["Input"];newlb=BasicAct[lb,action];
newstatecompl=CompleteState[newlb];
newstate=<|"Input"->newlb,"Terminal"->newstatecompl["Terminal"],"Value"->newstatecompl["Value"]|>;

(* compute reward *)
dvalue=newstatecompl["Value"]-state["Value"];
reward=Which[dvalue>0,dvalue,dvalue<=0,"RewardOffset"/.Options[LineBundleEnv]]+("StepPenalty"/.Options[LineBundleEnv]);

(* add boundary penalty *)
If[Max[newlb]== range[[2]] || Min[newlb]== range[[1]],reward=reward+("BoundaryPenalty"/.Options[LineBundleEnv])];
  
(* add terminal bonus *)
If[newstate["Terminal"],reward=reward+("TerminalBonus"/.Options[LineBundleEnv])];

<|"State"->state,"NewState"->newstate,"Action"->action,"Reward"->reward|>];


(* ::Section::Closed:: *)
(*End Package*)


End[]
EndPackage[]
