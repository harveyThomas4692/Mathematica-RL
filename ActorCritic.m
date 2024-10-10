(* ::Package:: *)

(* ::Title:: *)
(*ActorCritic: Policy-based RL algorithm*)


(* ::Section::Closed:: *)
(*Start package*)


BeginPackage["ActorCritic`"];


(* ::Section:: *)
(*Documentation*)


(* ::Subsection::Closed:: *)
(*General help*)


ActorCritic::usage="\[FilledSmallSquare]  This package implements an actor critic (A2C), policy-based RL algorithm.\n"<>
"\[FilledSmallSquare]  Global settings are contained in Options[ActorCritic] and can be changed with SetOptions.\n"<>
"\[FilledSmallSquare]  This package has to be combined with a realisation of an environment which provides modules to initialise a state and to act.\n"<>
"\[FilledSmallSquare]  The names of these modules are set in Options[ActorCritic], together with the dimensions of state tensor and action vector (one-hot encoded).\n"<>
"\[FilledSmallSquare]  States are associations of the form <|\"Input\" -> statetensor,\"Terminal\" -> True/False,\"Value\" -> value|>.\n"<>
"\[FilledSmallSquare]  The initialisation module has no mandatory input, so InitModuleName[options___], and its output has to be a state association.\n"<>
"\[FilledSmallSquare]  The act module should have the structure ActModuleName[stateassociation_,actionvector_,options___].\n"<>
"\[FilledSmallSquare]  The act module's output has to be of the form <|\"State\" -> state,\"NewState\"->newstate,\"Action\"->action,\"Reward\"->reward|> to represent a transition from a state to a newstate by an action, with the given reward.\n"<>
"\[FilledSmallSquare]  For a list of modules execute ?ActorCriticModules.";

ActorCriticModules::usage="Modules:\n"<>
"\[FilledSmallSquare]  Episode[initialstate,opt] creates an episode.\n"<>
"\[FilledSmallSquare]  PlotEpisode[episode,opt] plots rewards/values and state path of an episode.\n"<>
"\[FilledSmallSquare]  SampleEnv[batchsize,opt] to sample an enironment.\n"<>
"\[FilledSmallSquare]  ReLearn[opt] to carry out RL learning.";


(* ::Subsection::Closed:: *)
(*Module help*)


Status::usage="An association which containts the status during training, including the agents and a list of terminal states found during training."; 

Episode::usage="Episode[initialstate,opt] creates an episode.\n"<>
"\[FilledSmallSquare]  The initialstate should be given in the form <|\"Input\"->statevector,\"Terminal\"->True/False,\"Value\" -> value|>. Options:\n"<>
"\[FilledSmallSquare]  \"MaxEpisodeLength\" -> integer, specifies the maximum lenght of an episode (detault 32), realised unless a terminal state is encountered earlier.\n"<>
"\[FilledSmallSquare]  \"PolicyType\" -> string  is either \"probabilistic\" (default) or \"deterministic\" to determine actions using either the full probability vector or simply its maximal value.\n"<>
"\[FilledSmallSquare]  \"FullNet\"\[Rule]net, to provide a full AC network, as e.g. produced by ReLearn[] (default None).\n"<>
"\[FilledSmallSquare]  \"PolicyNet\"\[Rule]net, to provide a policy network (default None).\n"<>
"\[FilledSmallSquare]  \"ValueNet\"\[Rule]net, to provide a value network (default None).\n"<>
"\[FilledSmallSquare]  A policynet is required, provided either through the \"PolicyNet\" option or extracted from a full net provided with the \"FullNet\" option.\n"<>
"\[FilledSmallSquare]  \"Gamma\"\[Rule]real, to provide a discount factor (default 0.8).";

PlotEpisode::usage="PlotEpisode[episode] plots rewards/values and state path of an episode. Options:\n"<>
"\[FilledSmallSquare]  \"PlotLst\"->list, a list of integers specifying which quantities to plot (default {1,2,3} for reward, value and td return). Entries >3 indicate \"ValueLst\" quantities.\n"<>
"\[FilledSmallSquare]  \"Legend\"->list, to provide a legend for the \"ValueList\" quantities.\n"<>
"\[FilledSmallSquare]  \"ColorScheme\"->scheme, to specify the color scheme used for the plots.\n"<>
"\[FilledSmallSquare]  \"Labels\"->boolean, to switch labels on the state path on and off.";

SampleEnv::usage="SampleEnv[batchsize,opt] to sample an enironment. Options:\n"<>
"\[FilledSmallSquare]  \"EpisodeLength\"\[Rule]integer, to specify the maximal episode length for each agent (default 32). When this number or a terminal state is reached the agent is initialised.\n"<>
"\[FilledSmallSquare]  \"InitAgents\"\[Rule]boolean, to specify whether the complete status (including agents) should be initialised when module is called (default False).\n"<>
"\[FilledSmallSquare]  Networks should be provided using the same options as for Episode.";

ReLearn::usage="ReLearn[opt] to carry out RL learning.\n"<>
"\[FilledSmallSquare]  Trains policy network from environment. Output is an association which includes the NetTrain results object plus additional information. Options:\n"<>
"\[FilledSmallSquare]  \"FullNet\" -> net, to provide a fully initialised (or trained) AC net e.g. such as produced by ReLearn - this can be used to continue training a network (default None).\n"<>
"\[FilledSmallSquare]  \"PolicyNet\" -> list, to provide the rump of the policy net which will be completed with on output of the right size plus a softmax layer (default 3 FC layers, width 64).\n"<>
"\[FilledSmallSquare]  \"ValueNet\" -> list, to provide a rump value net which will be completed by adding a linear layer with a real output (default 3 FC layers, width 64).\n"<>
"\[FilledSmallSquare]  Unless a full net is provided by \"FullNet\", a complete network will be constructed and initialised from the policy and value nets provided by the previous two options.\n"<>
"\[FilledSmallSquare]  NetTrain options Method, MaxTrainingRound, Targetdevice and Batchsize.\n"<>
"\[FilledSmallSquare]  \"InitStatus\" -> boolean, if True (default) the variable Status will be initialised.\n"<>
"\[FilledSmallSquare]  \"LossScaling\"\[Rule]{p,v}, where p (v) scales the policyloss (valueloss) (default {0.5,1}).";


(* ::Section:: *)
(*Initialize*)


(* ::Subsection::Closed:: *)
(*Start-up messages*)


Begin["`Private`"]
Print[Style["ActorCritic: policy-based RL algorithm",Underlined,FontColor -> Blue,TextAlignment -> Center, FontSize -> 14,FontFamily->"Times"]];
Print[Style["Execute \!\(\*
StyleBox[\"?\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"ActorCritic\",\nFontWeight->\"Bold\"]\) for help",FontColor -> Black,TextAlignment -> Center, FontSize -> 14,FontFamily->"Times"]];


(* ::Subsection:: *)
(*Options*)


Options[ActorCritic]:={"InitModule"->LineBundleEnv`LineBundleInit,"ActModule"->LineBundleEnv`LineBundleAct,
                       "StateDim"->2,"ActionDim"->4,"NAgents"->4,"EnvironmentOptions"->LineBundleEnv`LineBundleEnv,
                       "ActionForm"->"Probabilities"};
Print[Style[Options[ActorCritic],"Output"]];  


(* ::Subsection::Closed:: *)
(*Global data*)


Status=<|"Agents"->None,"EpisodeLength"->{},"TerminalStates"->{},"NBatches"->0,"NEpisodes"->0,
         "TerminalList"->{},"FracTerminal"->{},"NTerminal"->{},"LastFracTerminal"->0,
         "EpisodeLengthList"->{},"AvEpisodeLength"->{},"LastAvEpisodeLength"->32|>;


(* ::Section:: *)
(*Modules*)


(* ::Subsection::Closed:: *)
(*Episode*)


Options[Episode]:={"MaxEpisodeLength"->32,"PolicyType"->"probabilistic","FullNet"->None,"PolicyNet"->None,
                   "ValueNet"->None,"Gamma"->0.8};

Episode[initialstate_,opt___]:=Module[
  {episodelength,policytype,actiondim,id,fullnet,policynet,valuenet,terminal,inputlst,valuelst,valuelstlst,
   actionlst,rewardlst,tdreturnlst,currentstate,n,actionprob,action,actionresult,gamma,tdreturn},

(* options *)
episodelength="MaxEpisodeLength"/.{opt}/.Options[Episode];
policytype="PolicyType"/.{opt}/.Options[Episode];
gamma="Gamma"/.{opt}/.Options[Episode];
actiondim="ActionDim"/.Options[ActorCritic];
id=IdentityMatrix[actiondim];

(* extract networks *)
fullnet="FullNet"/.{opt}/.Options[Episode];
If[fullnet===None,policynet="PolicyNet"/.{opt}/.Options[Episode]; valuenet="ValueNet"/.{opt}/.Options[Episode],
                  policynet=NetExtract[fullnet,"policynet"];valuenet=NetExtract[fullnet,"valuenet"]];
                  
(* initilise *)
terminal=initialstate["Terminal"]; inputlst={initialstate["Input"]}; valuelst={initialstate["Value"]};
If[KeyExistsQ[initialstate,"ValueLst"],valuelstlst={initialstate["ValueLst"]},valuelstlst={}];
actionlst={}; rewardlst={}; tdreturnlst={};
currentstate=initialstate; n=0;

(* main loop over steps of episode *)
While[!terminal && n<episodelength,

  (* determine action probability feeding current state into policy network *) 
  actionprob=policynet[currentstate["Input"]];

  (* determine action from probability vector *)
  If[policytype=="probabilistic",action=RandomSample[actionprob->id,1]//First,
                                 action=id[[Position[actionprob,Max[actionprob]][[1,1]]]]];
                                 
  (* act on current state to get next state and reward *)                            
  actionresult=("ActModule"/.Options[ActorCritic])[currentstate,action,opt];
  
  (* store new state, action and reward *)
  inputlst=Append[inputlst,actionresult["NewState"]["Input"]];
  actionlst=Append[actionlst,action];
  rewardlst=Append[rewardlst,actionresult["Reward"]];
  valuelst=Append[valuelst,actionresult["NewState"]["Value"]];
  If[KeyExistsQ[actionresult["NewState"],"ValueLst"],
     valuelstlst=Append[valuelstlst,actionresult["NewState"]["ValueLst"]]]; 
      
  (* if value network is present, compute td return *)
  If[valuenet=!=None,tdreturn=actionresult["Reward"]+gamma*valuenet[actionresult["NewState"]["Input"]];
                     tdreturnlst=Append[tdreturnlst,tdreturn]];  
  
  (* update *)
  terminal=actionresult["NewState"]["Terminal"];
  currentstate=actionresult["NewState"];
  n++];                                                                                                                                                                                                           
                                                                                                                                                                                                                                                                                                                                               
<|"Input"->inputlst,"Action"->actionlst,"Reward"->rewardlst,"TDReturn"->tdreturnlst,
  "Terminal"->terminal,"Value"->valuelst,"ValueLst"->valuelstlst|>];


(* ::Subsection::Closed:: *)
(*PlotEpisode*)


Options[PlotEpisode]:={"PlotLst"->{1,2,3},"ColorScheme"->"BrightBands","Legend"->{"1","2","3","4","5","6","7"},"Projection"->{},"Labels"->True};

PlotEpisode[episode_,opt___]:=Module[
{rewardlst,valuelst,statelst,i,ymin,ymax,xmax,rewardplot,statelen,v1,v2,statepts,stateplot,
 termtext,endpts,labels,xpts,ypts,xmin,xmargin,ymargin,valuelstlst,alllst,plotmarkers,plotstyle,plotlegends,
 plotlst,cs,sellst,legend,proj,plotlabels,tdreturnlst},

(* options *)
plotlst="PlotLst"/.{opt}/.Options[PlotEpisode];
cs=ColorData["ColorScheme"/.{opt}/.Options[PlotEpisode]];
legend="Legend"/.{opt}/.Options[PlotEpisode];
proj="Projection"/.{opt}/.Options[PlotEpisode];
plotlabels="Labels"/.{opt}/.Options[PlotEpisode];
rewardlst=episode["Reward"]; valuelst=episode["Value"]; statelst=episode["Input"];
tdreturnlst=episode["TDReturn"];
If[Length[episode["ValueLst"]]>0,valuelstlst=Transpose[episode["ValueLst"]],valuelstlst={}];
statelen=Length[Flatten[First[statelst]]];
If[proj=={},v1=Table[1,{statelen}]/Sqrt[statelen]; v2=Table[(-1)^i,{i,1,statelen}]/Sqrt[statelen],
           {v1,v2}=proj];

(* compute quantities for reward/value/return plot *)
If[Length[rewardlst]==0,rewardplot=None,
alllst=Join[{rewardlst,valuelst,tdreturnlst},valuelstlst];
sellst=alllst[[plotlst]];
ymin=Min[Flatten[sellst]]-1; ymax=Max[Flatten[sellst]]+1;
xmax=Max[Map[Length,sellst]];
If[episode["Terminal"],termtext=" (terminal)",termtext=" (non-terminal)"];
plotmarkers=Table[{"\[FilledCircle]",7},{Length[sellst]}];
plotstyle=Table[cs[i/Length[alllst]],{i,0,Length[alllst]}][[plotlst]];
plotlegends=Join[{"Reward","Value","TDReturn"},legend][[plotlst]];

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


Options[SampleEnv]:={"EpisodeLength"->32,"InitAgents"->False};

SampleEnv[batchsize_,opt___]:=Module[
  {episodelength,nagents,inputlst,actionlst,rewardlst,tdreturnlst,missing,eplength,agent,instate,episode,
   finalstate,eplenlst,agents,initagents,len},

(* options *)
episodelength="EpisodeLength"/.{opt}/.Options[SampleEnv];
nagents="NAgents"/.Options[ActorCritic];
initagents="InitAgents"/.{opt}/.Options[SampleEnv];

(* initialise agents if necessary *)
If[(Status["Agents"]==None) || Length[Status["Agents"]]!=nagents || initagents,
  Status=<|"Agents"->Table[("InitModule"/.Options[ActorCritic])[],{nagents}],
           "EpisodeLength"->Table[0,{nagents}],"TerminalStates"->{},"NBatches"->0,
           "NEpisodes"->0,"TerminalList"->{},"FracTerminal"->{},"NTerminal"->{},
           "LastFracTerminal"->0,"EpisodeLengthList"->{},"AvEpisodeLength"->{},
           "LastAvEpisodeLength"->episodelength|>];
  
(* check which agents have reached episode length and reset them *)
eplenlst=Status["EpisodeLength"]; agents=Status["Agents"];
For[agent=1,agent<=nagents,agent++,
  If[eplenlst[[agent]]>=episodelength,
    Status["EpisodeLengthList"]=Append[Status["EpisodeLengthList"],eplenlst[[agent]]];
    agents[[agent]]=("InitModule"/.Options[ActorCritic])[];
    eplenlst[[agent]]=0; Status["NEpisodes"]++;
    Status["TerminalList"]=Append[Status["TerminalList"],0]]];
Status["EpisodeLength"]=eplenlst; Status["Agents"]=agents;  
  
(* initialise lists *)
inputlst={}; actionlst={}; rewardlst={}; tdreturnlst={};

(* main loop over episodes *)
While[Length[inputlst]<batchsize,

  (* check how mach data is still missing and distribute onto agents *)
  missing=batchsize-Length[inputlst];
  eplength=Ceiling[missing/nagents];

  (* loop over agents *)
  For[agent=1,agent<=nagents,agent++,
  
    (* compute an episode *)
    instate=Status["Agents"][[agent]];
    episode=Episode[instate,"MaxEpisodeLength"->eplength,opt];
    
    (* update episode length *)
    eplenlst=Status["EpisodeLength"];
    eplenlst[[agent]]=eplenlst[[agent]]+Length[episode["Action"]];
    Status["EpisodeLength"]=eplenlst;
    
    (* add episode results to lists *)
    inputlst=Join[inputlst,Most[episode["Input"]]];
    actionlst=Join[actionlst,episode["Action"]];
    rewardlst=Join[rewardlst,episode["Reward"]];
    tdreturnlst=Join[tdreturnlst,episode["TDReturn"]];
    
    (* extract final state of episode *)
    finalstate=<|"Input"->Last[episode["Input"]],"Terminal"->episode["Terminal"],"Value"->Last[episode["Value"]]|>;
    
    (* update agent to random state if terminal and store terminal state *)
    eplenlst=Status["EpisodeLength"]; agents=Status["Agents"]; 
    If[finalstate["Terminal"],
      Status["EpisodeLengthList"]=Append[Status["EpisodeLengthList"],eplenlst[[agent]]];
      agents[[agent]]=("InitModule"/.Options[ActorCritic])[];
      Status["NEpisodes"]++; eplenlst[[agent]]=0;
      Status["TerminalList"]=Append[Status["TerminalList"],1];
      If[FreeQ[Status["TerminalStates"],finalstate],
        Status["TerminalStates"]=Append[Status["TerminalStates"],finalstate]],
      
      (* if state not terminal set agent to final state *)                                                    
      agents[[agent]]=finalstate];
    Status["EpisodeLength"]=eplenlst; Status["Agents"]=agents;
    
  ];
];      
  
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
   Status["AvEpisodeLength"]=Append[Status["AvEpisodeLength"],episodelength]];

(* store number of terminal states as a fct. of batch no. *)    
Status["NTerminal"]=Append[Status["NTerminal"],Length[Status["TerminalStates"]]]; 

(* increment batch counter *)
Status["NBatches"]++;

<|"Input"->Take[inputlst,batchsize],"Action"->Take[actionlst,batchsize],
  "TDReturn"->Take[tdreturnlst,batchsize]|>];         


(* ::Subsection:: *)
(*ReLearn*)


Options[ReLearn]:={Method->"ADAM",BatchSize->32,MaxTrainingRounds->1000,TargetDevice->"CPU","FullNet"->None,
                   "PolicyNet"->{64,ElementwiseLayer["SELU"],64,ElementwiseLayer["SELU"],64,ElementwiseLayer["SELU"]},
                   "ValueNet"->{64,ElementwiseLayer["SELU"],64,ElementwiseLayer["SELU"],64,ElementwiseLayer["SELU"]},
                   "InitStatus"->True,"LossScaling"->{0.5,1}};
                   
ReLearn[opt___]:=Module[
  {method,batchsize,maxtrainingrounds,targetdevice,fullnet,initstatus,statedim,actiondim,policylossnet,
   valuelossnet,policynet0,valuenet0,policynet,valuenet,tpm,netall,fracterminalplot,nterminalplot,fullres,
   lossscaling,avepisodeplot,xmax,ymax,actionform},

(* options *)
method=Method/.{opt}/.Options[ReLearn];
batchsize=BatchSize/.{opt}/.Options[ReLearn];
maxtrainingrounds=MaxTrainingRounds/.{opt}/.Options[ReLearn];
targetdevice=TargetDevice/.{opt}/.Options[ReLearn];
fullnet="FullNet"/.{opt}/.Options[ReLearn];
policynet0="PolicyNet"/.{opt}/.Options[ReLearn];
valuenet0="ValueNet"/.{opt}/.Options[ReLearn];
initstatus="InitStatus"/.{opt}/.Options[ReLearn];
lossscaling="LossScaling"/.{opt}/.Options[ReLearn];
actionform="ActionForm"/.{opt}/.Options[ActorCritic];

(* initialise status if required *)
If[initstatus,Status=<|"Agents"->None,"EpisodeLength"->{},"TerminalStates"->{},"NBatches"->0,"NEpisodes"->0, 
                       "TerminalList"->{},"FracTerminal"->{},"NTerminal"->{},"LastFracTerminal"->0,
                       "EpisodeLengthList"->{},"AvEpisodeLength"->{},"LastAvEpisodeLength"->batchsize|>];

(* build a new network from given rump policy and value networks and initialise *)                       
If[fullnet===None,

  (* dimensions *)
  statedim="StateDim"/.Options[ActorCritic];
  actiondim="ActionDim"/.Options[ActorCritic];
  
  (* complete policy net *) 
  policynet= NetChain[Join[policynet0,{actiondim,SoftmaxLayer[]}],"Input"->statedim];
  
  (* complete value net *)  
  valuenet=NetChain[Join[valuenet0,{{}}],"Input"->statedim]; 
  
  (* define policy lossnet *) 
  policylossnet=NetGraph[{CrossEntropyLossLayer[actionform],ThreadingLayer[#1*#2&]},
                 {NetPort["TDReturn"]->2,{NetPort["Input"],NetPort["Action"]}->1->2->NetPort["PolicyLoss"]}];
   
  (* define value lossnet *)
  valuelossnet=NetGraph[{MeanSquaredLossLayer[]},{NetPort["Input"]->1->NetPort["ValueLoss"],NetPort["TDReturn"]->1}];  
 
  (* combine to complete network *)
  fullnet=NetGraph[<|"policynet"->policynet,"valuenet"->valuenet,
                     "policylossnet"->policylossnet,"valuelossnet"->valuelossnet|>,
                     {NetPort["policynet","Output"]->NetPort["policylossnet","Input"],
                      NetPort["valuenet","Output"]->NetPort["valuelossnet","Input"]}]//NetInitialize;
 
];

(* train network *)
tpm=<|"Measurement"->NetPort["TDReturn"],"Aggregation"->"Mean"|>;
netall=Monitor[

         NetTrain[fullnet,SampleEnv[#BatchSize,"FullNet"->#Net,opt]&,All,
                  LossFunction->{"PolicyLoss"->Scaled[lossscaling[[1]]],"ValueLoss"->Scaled[lossscaling[[2]]]},
                  Method->method,MaxTrainingRounds->maxtrainingrounds,BatchSize->batchsize,TargetDevice->targetdevice,
                  TrainingProgressMeasurements->tpm],
                  
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

(* association containing complete result *)
fullres=<|"NetAll"->netall,
          "Method"->method,
          "EnvironmentOptions"->Options[("EnvironmentOptions"/.Options[ActorCritic])],
          "InitState"->("InitState"/.{opt}/.Options[SampleEnv]),
          "MaxEpisodeLength"->("MaxEpisodeLength"/.{opt}/.Options[Episode]),
          "PolicyType"->("PolicyType"/.{opt}/.Options[Episode]),
          "Gamma"->("Gamma"/.{opt}/.Options[SampleEnv]),
          "FullNet"->netall["TrainedNet"],
          "PolicyNet"->NetExtract[netall["TrainedNet"],"policynet"],
          "ValueNet"->NetExtract[netall["TrainedNet"],"valuenet"],
          "Loss"->netall["FinalPlots"]["Loss"],
          "PolicyLoss"->netall["FinalPlots"]["PolicyLoss"],
          "ValueLoss"->netall["FinalPlots"]["ValueLoss"],
          "TDReturn"->netall["FinalPlots"]["TDReturn"],
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
