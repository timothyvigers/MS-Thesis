library(DiagrammeR)
med1 = grViz("
  digraph{
  node [shape = box,fontname = Helvetica]
  b[label= 'Metabolite']
  a[label= 'Methylation']
  c[label= 'IA']
  a -> b 
  b -> c 
  {rank = same;a -> c}
  }
")
med1

med2 = grViz("
  digraph{
  node [shape = box,fontname = Helvetica]
  b[label= 'Methylation']
  a[label= 'Metabolite']
  c[label= 'IA']
  a -> b 
  b -> c 
  {rank = same;a -> c}
  }
")
med2

flow = grViz("digraph {

graph [layout = dot, rankdir = LR]

node [shape = rectangle]

birth [label = 'Birth']
first [label = 'First Abx+ \n Visit']
second [label = 'Second Abx+ \n Visit (IA)']
psv [label = 'PSV']
sv [label = 'SV']
t1d [label = 'T1D']

age9 [label =  '9 Month Visit']
age15 [label =  '15 Month Visit']
age24 [label =  '24 Month Visit']
annual [label =  'Annual Visits']
accelerated [label =  'Accelerated Protocol']

# edge definitions with the node IDs
birth -> first -> second -> t1d;
birth -> age9 -> age15 -> age24 -> annual;
{rank=same; second -> sv}
{rank=same; first -> psv}
sv -> accelerated
}")
flow


med = grViz("
  digraph{
  node [shape = plain]
  M
  A
  Y
  A -> M 
  M -> Y 
  {rank = same;A -> Y}
  }
")
med
