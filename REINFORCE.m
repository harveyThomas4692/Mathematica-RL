(* ::Package:: *)

(* ::Title:: *)
(*REINFORCE: Policy-based RL algorithm*)


(* ::Section::Closed:: *)
(*Start package*)


BeginPackage["REINFORCE`"];


(* ::Section:: *)
(*Documentation*)


(* ::Subsection::Closed:: *)
(*General help*)


REINFORCE::usage="\[FilledSmallSquare]  This package implements the REINFORCE, policy-based RL algorithm.\n"<>
"\[FilledSmallSquare]  Global settings are contained in Options[REINFORCE] and can be changed with SetOptions.\n"<>
"\[FilledSmallSquare]  This package has to be combined with a realisation of an environment which provides modules to initialise a state and to act.\n"<>
"\[FilledSmallSquare]  The names of these modules are set in Options[REINFORCE], together with the dimensions of state tensor and action.\n"<>
"\[FilledSmallSquare]  The action shouould be provided as a one-hot vector (for \"ActionForm\"->\"Probabilities\") or as an integer (for \"ActionForm\"->\"Index\").\n"<>
"\[FilledSmallSquare]  States are associations of the form <|\"Input\" -> statetensor,\"Terminal\" -> True/False,\"Value\" -> value|>.\n"<>
"\[FilledSmallSquare]  The initialisation module has no mandatory input, so InitModuleName[options___], and its output has to be a state association.\n"<>
"\[FilledSmallSquare]  The act module should have the structure ActModuleName[stateassociation_,actionvector_,options___].\n"<>
"\[FilledSmallSquare]  The act module's output has to be of the form <|\"State\" -> state,\"NewState\"->newstate,\"Action\"->action,\"Reward\"->reward|> to represent a transition from a state to a newstate by an action, with the given reward.\n"<>
"\[FilledSmallSquare]  For a list of modules execute ?REINFORCEModules.";

REINFORCEModules::usage="Modules:\n"<>
"\[FilledSmallSquare]  Episode[initialstate,policynet,opt] creates an episode.\n"<>
"\[FilledSmallSquare]  PlotEpisode[episode] plots rewards/values and state path of an episode.\n"<>
"\[FilledSmallSquare]  SampleEnv[policynet,batchsize,opt] to sample an enironment.\n"<>
"\[FilledSmallSquare]  ReLearn[opt] to carry out RL learning.";


(* ::Subsection::Closed:: *)
(*Module help*)


Status::usage="An association which containts the status during training, including the agent and a list of terminal states found during training."; 

Episode::usage="Episode[initialstate,policynet,opt] creates an episode.\n"<>
"\[FilledSmallSquare]  The initialstate should be given in the form <|\"Input\"->statevector,\"Terminal\"->True/False,\"Value\" -> value|>. Options:\n"<>
"\[FilledSmallSquare]  \"MaxEpisodeLength\" -> integer, specifies the maximum lenght of an episode (detault 16), realised unless a terminal state is encountered earlier.\n"<>
"\[FilledSmallSquare]  \"PolicyType\" -> string  is either \"probabilistic\" (default) or \"deterministic\" to determine actions using either the full probability vector or simply its maximal value.\n";

PlotEpisode::usage="PlotEpisode[episode,opt] plots rewards/values and state path of an episode. Options:\n"<>
"\[FilledSmallSquare]  \"PlotLst\"->list, a list of integers specifying which quantities to plot (default {1,2} for reward and value). Entries >2 indicate \"ValueLst\" quantities.\n"<>
"\[FilledSmallSquare]  \"Legend\"->list, to provide a legend for the \"ValueList\" quantities.\n"<>
"\[FilledSmallSquare]  \"ColorScheme\"->scheme, to specify the color scheme used for the plots.\n"<>
"\[FilledSmallSquare]  \"Labels\"->boolean, to switch labels on the state path on and off.";

SampleEnv::usage="SampleEnv[net,batchsize,opt] to sample an enironment. Options:\n"<>
"\[FilledSmallSquare]  \"InitState\" -> boolean, to specify if each episode should start from a new random state (default False).\n"<>
"\[FilledSmallSquare]  \"Gamma\" -> real, specifies the discount rate (default 0.98).\n"<>
"\[FilledSmallSquare]  \"PolicyNet\" -> boolean, True (default) if net is policynet, False if policynet is part of a larger net.";

ReLearn::usage="ReLearn[opt] to carry out RL learning.\n"<>
"\[FilledSmallSquare]  Trains policy network from environment. Output is an association which includes the NetTrain results object plus additional information. Options:\n"<>
"\[FilledSmallSquare]  NetTrain options Method, MaxTrainingRound, Targetdevice and Batchsize.\n"<>
"\[FilledSmallSquare]  \"InputNet\" -> net, to specify the rump policynet, with final layer and softmax added automatically if \"InitNet\" -> True (default 3 FC layers, width 64).\n"<>
"\[FilledSmallSquare]  \"InitNet\" -> boolean, if True (default) network is contructed from InputNet and initialised, otherwise it is assumed InputNet is a fully initialised (or trained) RL network.\n"<>
"\[FilledSmallSquare]  \"InitStatus\" -> boolean, if True (default) the variable Status will be initialised.";


(* ::Section:: *)
(*Initialize*)


(* ::Subsection::Closed:: *)
(*Start-up messages*)


Begin["`Private`"]
Print[Style["REINFORCE: Policy-based RL algorithm",Underlined,FontColor -> Blue,TextAlignment -> Center, FontSize -> 14,FontFamily->"Times"]];
Print[Style["Execute \!\(\*
StyleBox[\"?\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"REINFORCE\",\nFontWeight->\"Bold\"]\) for help",FontColor -> Black,TextAlignment -> Center, FontSize -> 14,FontFamily->"Times"]];


(* ::Subsection::Closed:: *)
(*Options*)


Options[REINFORCE]:={"InitModule"->FNEnv`FNInit,"ActModule"->FNEnv`FNAct,"StateDim"->10,"ActionDim"->20,
                     "EnvironmentOptions"->FNEnv`FNEnv,"ActionForm"->"Probabilities"};
Print[Style[Options[REINFORCE],"Output"]];


(* ::Subsection::Closed:: *)
(*Global data*)


(* REINFORCE`Agent=None;
REINFORCE`TerminalStates={}; *)
Status=<|"Agent"->None,"TerminalStates"->{},"NEpisodes"->0,"TerminalList"->{},
         "FracTerminal"->{0},"NTerminal"->{},"LastFracTerminal"->0,
         "EpisodeLengthList"->{},"AvEpisodeLength"->{},"LastAvEpisodeLength"->32,"LastBatch"->None|>;


(* ::Section:: *)
(*Modules*)


(* ::Subsection::Closed:: *)
(*Episode*)


Options[Episode]:={"MaxEpisodeLength"->32,"PolicyType"->"probabilistic"};

Episode[initialstate_,policynet_,opt___]:=Module[
  {episodelength,n,inputlst,actionlst,rewardlst,terminal,action,actionresult,currentstate,actiondim,id,policytype,
   actionprob,valuelst,valuelstlst},

(* options *)
episodelength="MaxEpisodeLength"/.{opt}/.Options[Episode];
policytype="PolicyType"/.{opt}/.Options[Episode];
actiondim="ActionDim"/.Options[REINFORCE];
id=IdentityMatrix[actiondim];

(* initilise *)
terminal=initialstate["Terminal"]; inputlst={initialstate["Input"]}; actionlst={}; rewardlst={};
valuelst={initialstate["Value"]};
If[KeyExistsQ[initialstate,"ValueLst"],valuelstlst={initialstate["ValueLst"]},valuelstlst={}];
currentstate=initialstate; n=0;

(* main loop over steps of episode *)
While[!terminal && n<episodelength,

  (* determine action probability feeding current state into policy network *) 
  actionprob=policynet[Last[inputlst]];
  
  (* determine action from probability vector *)
  If[policytype=="probabilistic",action=RandomSample[actionprob->id,1]//First,
                                 action=id[[Position[actionprob,Max[actionprob]][[1,1]]]]];
                                 
  (* act on current state to get next state and reward *)                            
  actionresult=("ActModule"/.Options[REINFORCE])[currentstate,action,opt];
  
  (* store new state, action and reward *)
  inputlst=Append[inputlst,actionresult["NewState"]["Input"]];
  actionlst=Append[actionlst,action];
  rewardlst=Append[rewardlst,actionresult["Reward"]];
  valuelst=Append[valuelst,actionresult["NewState"]["Value"]];
  If[KeyExistsQ[actionresult["NewState"],"ValueLst"],
     valuelstlst=Append[valuelstlst,actionresult["NewState"]["ValueLst"]]];
  
  (* update *)
  terminal=actionresult["NewState"]["Terminal"];
  currentstate=actionresult["NewState"];
  n++];
  
<|"Input"->inputlst,"Action"->actionlst,"Reward"->rewardlst,"Terminal"->terminal,
  "Value"->valuelst,"ValueLst"->valuelstlst|>];


(* ::Subsection::Closed:: *)
(*PlotEpisode*)


Options[PlotEpisode]:={"PlotLst"->{1,2},"ColorScheme"->"BrightBands","Legend"->{"1","2","3","4","5","6","7"},"Projection"->{},"Labels"->True};

PlotEpisode[episode_,opt___]:=Module[
{rewardlst,valuelst,statelst,i,ymin,ymax,xmax,rewardplot,statelen,v1,v2,statepts,stateplot,
 termtext,endpts,labels,xpts,ypts,xmin,xmargin,ymargin,valuelstlst,alllst,plotmarkers,plotstyle,plotlegends,
 plotlst,cs,sellst,legend,proj,plotlabels},

(* options *)
plotlst="PlotLst"/.{opt}/.Options[PlotEpisode];
cs=ColorData["ColorScheme"/.{opt}/.Options[PlotEpisode]];
legend="Legend"/.{opt}/.Options[PlotEpisode];
proj="Projection"/.{opt}/.Options[PlotEpisode];
plotlabels="Labels"/.{opt}/.Options[PlotEpisode];
rewardlst=episode["Reward"]; valuelst=episode["Value"]; statelst=episode["Input"];
If[Length[episode["ValueLst"]]>0,valuelstlst=Transpose[episode["ValueLst"]],valuelstlst={}];
statelen=Length[Flatten[First[statelst]]];
If[proj=={},v1=Table[1,{statelen}]/Sqrt[statelen]; v2=Table[(-1)^i,{i,1,statelen}]/Sqrt[statelen],
           {v1,v2}=proj];

(* compute quantities for reward/value plot *)
If[episode["Terminal"],termtext=" (terminal)",termtext=" (non-terminal)"];
If[Length[rewardlst]==0,rewardplot=None,
alllst=Join[{rewardlst,valuelst},valuelstlst];
sellst=alllst[[plotlst]];
ymin=Min[Flatten[sellst]]-1; ymax=Max[Flatten[sellst]]+1;
xmax=Max[Map[Length,sellst]];
plotmarkers=Table[{"\[FilledCircle]",7},{Length[sellst]}];
plotstyle=Table[cs[i/Length[alllst]],{i,0,Length[alllst]}][[plotlst]];
plotlegends=Join[{"Reward","Value"},legend][[plotlst]];

(* create reward/value plot *)
rewardplot=ListLinePlot[sellst,Frame->True,PlotRange->{{1,xmax},{ymin,ymax}},
             PlotMarkers->plotmarkers,PlotStyle->plotstyle,PlotLegends->plotlegends,GridLines->Automatic,
             GridLinesStyle->Directive[Thin,Gray],PlotLabel->"episode characteristics"<>termtext,
             ImageSize->Medium];
];
                          
(* compute quantities for state plot *)           
statepts=Table[{v1.Flatten[statelst[[i]]],v2.Flatten[statelst[[i]]]},{i,1,Length[statelst]}]//N;
endpts={{First[statepts]},{Last[statepts]}};
If[plotlabels,labels=Round[episode["Value"],0.1],labels={}];
xpts=Map[First,statepts]; ypts=Map[Last,statepts];
xmin=Min[xpts];xmax=Max[xpts];ymin=Min[ypts];ymax=Max[ypts];
xmargin=(xmax-xmin)/10;ymargin=(ymax-ymin)/10;
xmin=xmin-xmargin;xmax=xmax+xmargin;ymin=ymin-ymargin;ymax=ymax+ymargin;

stateplot=Show[ListLinePlot[statepts->labels,Frame->True,PlotStyle->cs[0.5],PlotMarkers->{"\[FilledCircle]",7},GridLines->Automatic,
               GridLinesStyle->Directive[Thin,Gray],PlotLabel->"2d projection of state path"<>termtext,ImageSize->Medium,
               PlotRange->{{xmin,xmax},{ymin,ymax}}],
               ListPlot[endpts,PlotStyle->{{cs[1],PointSize[0.03]},{cs[0],PointSize[0.03]}},PlotLegends->{"Start","End"}]];             
                          
<|"Rewards"->rewardplot,"States"->stateplot|>];


(* ::Subsection::Closed:: *)
(*SampleEnv*)


Options[SampleEnv]:={"InitState"->True,"Gamma"->0.98,"PolicyNet"->True};

SampleEnv[net_,batchsize_,opt___]:=Module[
  {initialstate,inputlst,actionlst,rewardlst,episode,initstate,gamma,randomstate,returnlst,
   finalstate,policynet,netstruct,termlen,eplength,len,batch},

(* options *)
initstate="InitState"/.{opt}/.Options[SampleEnv];
netstruct="PolicyNet"/.{opt}/.Options[SampleEnv];
gamma="Gamma"/.{opt}/.Options[SampleEnv];
If[netstruct,policynet=net,policynet=NetExtract[net,"policynet"]];

(* determine initial state *)
If[(Status["Agent"]==None) || initstate,
  Status["Agent"]=("InitModule"/.Options[REINFORCE])[]];

(* initialise lists *)
inputlst={}; actionlst={}; rewardlst={}; returnlst={}; 

(* main loop over episodes *)
While[Length[inputlst]<batchsize,

  initialstate=Status["Agent"];
   
  (* carry out episode *) 
  episode=Episode[initialstate,policynet,opt,"MaxEpisodeLength"->(batchsize)];
  Status["NEpisodes"]++;
  
  (* add episode data to lists *)
  inputlst=Join[inputlst,Most[episode["Input"]]];
  actionlst=Join[actionlst,episode["Action"]];
  rewardlst=Join[rewardlst,episode["Reward"]];
  eplength=Length[episode["Action"]];
  
  (* extract final state of episode *)
  finalstate=<|"Input"->Last[episode["Input"]],"Terminal"->episode["Terminal"],"Value"->Last[episode["Value"]]|>;
  
  (* add episode length to status *)
  Status["EpisodeLengthList"]=Append[Status["EpisodeLengthList"],eplength];
  
  (* if episode is terminal... *)
  If[episode["Terminal"],
    Status["TerminalList"]=Append[Status["TerminalList"],1];
    Status["Agent"]=("InitModule"/.Options[REINFORCE])[];
    If[FreeQ[Status["TerminalStates"],finalstate],Status["TerminalStates"]=Append[Status["TerminalStates"],finalstate]];
    returnlst=Join[returnlst,REINFORCE`Private`ComputeReturns[Take[rewardlst,{Length[returnlst]+1,Length[rewardlst]}],gamma]],
    
  (* if episode is not terminal *)
    Status["TerminalList"]=Append[Status["TerminalList"],0];
    Status["Agent"]=finalstate]  
                      
];

(* calculate returns which haven't yet been calculated *)
returnlst=Join[returnlst,REINFORCE`Private`ComputeReturns[Take[rewardlst,{Length[returnlst]+1,Length[rewardlst]}],gamma]];

(* maintain terminal list: only keep the last 100 episodes, store fraction of terminal episodes as a fct. of batch no. *)
If[Length[Status["TerminalList"]]>=1,
  len=Min[Length[Status["TerminalList"]],100];  
  Status["TerminalList"]=Take[Status["TerminalList"],-len];
  Status["FracTerminal"]=Append[Status["FracTerminal"],Mean[Status["TerminalList"]]//N];
  Status["LastFracTerminal"]=Round[Last[Status["FracTerminal"]],0.01],
  Status["FracTerminal"]=Append[Status["FracTerminal"],0]];

(* maintain episode length list in a similar way *)
If[Length[Status["EpisodeLengthList"]]>=1,
   len=Min[Length[Status["EpisodeLengthList"]],100];
   Status["EpisodeLengthList"]=Take[Status["EpisodeLengthList"],-len];
   Status["AvEpisodeLength"]=Append[Status["AvEpisodeLength"],Mean[Status["EpisodeLengthList"]]//N];
   Status["LastAvEpisodeLength"]=Round[Last[Status["AvEpisodeLength"]],0.1],
   Status["AvEpisodeLength"]=Append[Status["AvEpisodeLength"],batchsize]];

(* store number of terminal states as a fct. of batch no. *)    
Status["NTerminal"]=Append[Status["NTerminal"],Length[Status["TerminalStates"]]]; 

batch=<|"Input"->Take[inputlst,batchsize],"Action"->Take[actionlst,batchsize],"Return"->Take[returnlst,batchsize]|>;
Status["LastBatch"]=batch;

batch];


(* module to compute returns from rewards *)
REINFORCE`Private`ComputeReturns[rewards_,gamma_]:=Module[{returns,i,n,returnsnorm},

n=Length[rewards];
If[n==0,Return[{}]];
returns={Last[rewards]};
For[i=n-1,i>0,i--,returns=Prepend[returns,rewards[[i]]+gamma*First[returns]]];
returnsnorm=Table[returns[[i]]/(1-gamma^(n-i+1)),{i,1,n}];

returnsnorm];


(* ::Subsection::Closed:: *)
(*ReLearn*)


Options[ReLearn]:={Method->"ADAM",BatchSize->32,MaxTrainingRounds->1000,TargetDevice->"CPU",
                   "InputNet"->{64,ElementwiseLayer["SELU"],64,ElementwiseLayer["SELU"],64,ElementwiseLayer["SELU"]},
                   "InitNet"->True,"InitStatus"->True};

ReLearn[opt___]:=Module[
  {method,batchsize,maxtrainingrounds,inputnet,targetdevice,netname,statedim,actiondim,lossnet,policynet,
   fullnet,tpm,netall,init,initstatus,measurenet,fracterminalplot,fullres,nterminalplot,xmax,ymax,
   avepisodeplot,actionform},

(* options *)
method=Method/.{opt}/.Options[ReLearn];
batchsize=BatchSize/.{opt}/.Options[ReLearn];
maxtrainingrounds=MaxTrainingRounds/.{opt}/.Options[ReLearn];
targetdevice=TargetDevice/.{opt}/.Options[ReLearn];
inputnet="InputNet"/.{opt}/.Options[ReLearn];
init="InitNet"/.{opt}/.Options[ReLearn];
initstatus="InitStatus"/.{opt}/.Options[ReLearn];
actionform="ActionForm"/.{opt}/.Options[REINFORCE];

(* initialise status if required *)
If[initstatus,Status=<|"Agent"->None,"TerminalStates"->{},"NEpisodes"->0,"TerminalList"->{},
         "FracTerminal"->{0},"NTerminal"->{},"LastFracTerminal"->0,
         "EpisodeLengthList"->{},"AvEpisodeLength"->{},"LastAvEpisodeLength"->batchsize,"LastBatch"->None|>];

(* build a new network from a given rump policy network and initialise *)
If[init,

(* dimensions *)
statedim="StateDim"/.Options[REINFORCE];
actiondim="ActionDim"/.Options[REINFORCE];

(* define lossnet *) 
lossnet=NetGraph[{CrossEntropyLossLayer[actionform],ThreadingLayer[#1*#2&]},
                 {NetPort["Return"]->2,{NetPort["Input"],NetPort["Action"]}->1->2->NetPort["Loss"]}];
                 
(* complete policy net *) 
policynet= NetChain[Join[inputnet,{actiondim,SoftmaxLayer[]}],"Input"->statedim];                  

(* combine policy net and loss net *)
fullnet=NetGraph[<|"policynet"->policynet,"lossnet"->lossnet|>,{NetPort["policynet","Output"]->NetPort["lossnet","Input"]}]//NetInitialize,

(* or, if a full, initialised network is provided, use this *)

fullnet=inputnet];

(* train network *)
tpm=<|"Measurement"->NetPort["Return"],"Aggregation"->"Mean"|>;
netall=Monitor[

         NetTrain[fullnet,SampleEnv[#Net,#BatchSize,"PolicyNet"->False,opt]&,All,LossFunction->"Loss",
                  Method->method,MaxTrainingRounds->maxtrainingrounds,BatchSize->batchsize,
                  TargetDevice->targetdevice,TrainingProgressMeasurements->tpm],
        
       Style["# terminal states: "<>ToString[Length[Status["TerminalStates"]]]<>
       ",   terminal episodes: "<>ToString[Status["LastFracTerminal"]]<>
       ",   episode length: "<>ToString[Status["LastAvEpisodeLength"]],FontColor -> Black, FontSize -> 12,FontFamily->"Times"]];

(* final plots *)
xmax=Length[Status["FracTerminal"]]+1;
fracterminalplot=ListLinePlot[Status["FracTerminal"],Frame->True,PlotStyle->Orange,GridLines->Automatic,
                          ImageSize->Scaled[0.3],FrameLabel->{"rounds","terminal fraction"},PlotRange->{{0,xmax},{0,1}}];
xmax=Length[Status["NTerminal"]]+1; ymax=Max[Status["NTerminal"]]+1;
nterminalplot=ListLinePlot[Status["NTerminal"],Frame->True,PlotStyle->Orange,GridLines->Automatic,
                          ImageSize->Scaled[0.3],FrameLabel->{"rounds","# terminal states"},PlotRange->{{0,xmax},{0,ymax}}];
xmax=Length[Status["AvEpisodeLength"]]+1;ymax=Max[Status["AvEpisodeLength"]]+1;                         
avepisodeplot=ListLinePlot[Status["AvEpisodeLength"],Frame->True,PlotStyle->Orange,GridLines->Automatic,
                          ImageSize->Scaled[0.3],FrameLabel->{"rounds","episode length"},PlotRange->{{0,xmax},{0,ymax}}];

fullres=<|"NetAll"->netall,
          "Method"->method,
          "EnvironmentOptions"->Options[("EnvironmentOptions"/.Options[REINFORCE])],
          "InitState"->("InitState"/.{opt}/.Options[SampleEnv]),
          "MaxEpisodeLength"->("MaxEpisodeLength"/.{opt}/.Options[Episode]),
          "PolicyType"->("PolicyType"/.{opt}/.Options[Episode]),
          "Gamma"->("Gamma"/.{opt}/.Options[SampleEnv]),
          "TrainedNet"->netall["TrainedNet"],
          "PolicyNet"->NetExtract[netall["TrainedNet"],"policynet"],
          "Loss"->netall["FinalPlots"]["Loss"],
          "Return"->netall["FinalPlots"]["Return"],
          "FracTerminal"->fracterminalplot,
          "NTerminal"->nterminalplot,
          "AvEpisodeLength"->avepisodeplot,
          "NTerminalStates"->Length[Status["TerminalStates"]],
          "TerminalStates"->Status["TerminalStates"]|>;


fullres];


(* ::Section::Closed:: *)
(*End Package*)


End[]
EndPackage[]
