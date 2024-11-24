import math
import matplotlib.pyplot as plt
import numpy as np
from typing import List,Callable

epsilon=1e-3
memo={}

#Point class 정의
class Point:
    def __init__(self, coordinates:List[float]):
        self.coordinates=coordinates

    def add(self,other:'Point')->'Point':
        dim=len(self.coordinates)
        result=[self.coordinates[i]+other.coordinates[i] for i in range(dim)]
        return Point(result)
    
    def scale(self, k:float)->'Point':
        result=[k*x for x in self.coordinates]
        return Point(result)
    

    def inner(self, other: 'Point')->float:
        dim=len(self.coordinates)
        return sum(self.coordinates[i]*other.coordinates[i] for i in range(dim))
    
    def outer(self,other:'Point')->'Point': # 3차원 벡터의 외적
        if self.dimention()!=3 and other.dimention()!=3:
            raise ValueError("outer must can be in 3 dimention")
        x1,y1,z1=self.coordinates()
        x2,y2,z2=other.coordinates()
        result=[y1*z2-z1*y2,z1*x2-x1*z2,x1*y2-y1*x2]
        return Point(result)

    def norm(self)->float:
        result=0
        for i in range(self.dimention()):
            result+=self.coordinates[i]**2
        return math.sqrt(result)
    
    def degree(self,other:'Point')->float:
        costheta=self.inner(other)/(self.norm()*other.norm())
        theta=math.acos(costheta)
        return theta
    
    def distance(self, other: 'Point') -> float:
        diff = Point([self.coordinates[i] - other.coordinates[i] for i in range(len(self.coordinates))])
        return math.sqrt(diff.inner(diff))
    
    def dimention(self):
        return len(self.coordinates)
    
    def __repr__(self):
        return f"point({self.coordinates})"
    


# Set, FiniteSet, InfiniteSet class 정의
class Set():
    def __init__(self,element:list):
        self.element=element

class FiniteSet(Set):
    def __init__(self,element:list):
        super().__init__(element)
        self.element=element

    def len(self,element:list):
        return len(element)
    
    def Contain(self,urea):
        if urea in self.element:
            return True
        else:
            return False
    
class InfiniteSet():
    def __init__(self,lowerbound:Point,upperbound:Point,condition:Callable[[Point],bool]):
        if lowerbound.dimention()!=upperbound.dimention():
            raise ValueError("The dimensions of the boundary values ​​must be the same")
        self.lowerbound=lowerbound
        self.upperbound=upperbound
        self.condition=condition

    def Contain(self,point)->bool:
        return self.condition(self.lowerbound,self.upperbound,point)

   


#Domain class 정의
class Domain(InfiniteSet):
    def __init__(self,lowerbound:Point,upperbound:Point,condition:Callable[[Point],bool]):
        super().__init__(lowerbound,upperbound,condition)
        if lowerbound.dimention()!=upperbound.dimention():
            raise ValueError("The dimensions of the boundary values ​​must be the same")
        self.lowerbound=lowerbound
        self.upperbound=upperbound
        self.condition=condition

    def dimention(self) -> int:
        return self.lowerbound.dimention()
        

        
    def Bounds(self)->List[Point]:
        return [self.lowerbound,self.upperbound]
    
    def __repr__(self):
        return f"Domain(lowerbound={self.lowerbound},upperbound={self.upperbound})"
    

#Function class 정의
class Function(Domain):
    def __init__(self,lowerbound:Point,upperbound:Point,mapping:Callable[[Point],float],condition:Callable[[Point],bool]):
        super().__init__(lowerbound,upperbound,condition)
        self.mapping=mapping
        self.lowerbound=lowerbound
        self.upperbound=upperbound
        self.condition=condition

    def Map(self,point:Point)->float:
        if point.dimention()!=self.dimention():
            raise ValueError("Point dimention does not match function dimention")
        value=self.mapping(point)
        result=Point(point.coordinates)
        if self.Contain(result):
            return value
        else:
            raise ValueError("Point is not contained in domain")
        


#Parameter class 정의
class Parameter(Domain):
    def __init__(self,lowerbound:Point,upperbound:Point,mapping:Callable[[Point],Point]):
        super().__init__(lowerbound,upperbound,condition)
        self.mapping=mapping
        self.lowerbound=lowerbound
        self.upperbound=upperbound
        self.condition=condition

    def Map(self,point:Point)->'Point':
        if point.dimention()!=self.dimention():
            raise ValueError("Point dimention does not match parameter dimention")
        result=self.mapping(point)
        return result
    

#간단한 diffequal(미분방정식) class 정의
class diffequal():
    def __init__(self,lowerbound:Point,upperbound:Point):
        self.lowerbound=lowerbound
        self.upperbound=upperbound
    def Map(t,y):
        return -y
    


# def diff(편미분) 
def diff(func:Function,point:Point,Parameter)->float:
    if Parameter>func.dimention():
        raise ValueError("Dimention is not matched")
    realfuncvalue=func.Map(point)
    point.coordinates[Parameter]+=epsilon
    right_funcvalue=func.Map(point)
    point.coordinates[Parameter]-=2*epsilon
    left_funcvalue=func.Map(point)
    if abs(right_funcvalue-left_funcvalue)>epsilon*1000:
        raise ValueError("Non-continous in point")
    right_diffvalue=(right_funcvalue-realfuncvalue)/epsilon
    left_diffvalue=(realfuncvalue-left_funcvalue)/epsilon
    if abs(right_diffvalue-left_diffvalue)<epsilon*1000:
        return (right_diffvalue+left_diffvalue)/2
    else:
        raise ValueError("Non-differentiable in point")
    

# def gradient
def gradient(func:Function,point:Point)->'Point':
    result=[]
    for i in range(func.dimention()):
        result.append(diff(func,point,i))
    grt=Point(result)
    return grt


#def newton_labson(뉴천-랩슨법, 방정식의 근 구하기)
def newton_labson(func:Function,point:Point,number:int):
    if func.dimention()!=1 or point.dimention()!=1:
        raise ValueError("newton-labson must can be in 2 dimention")
    def newton(n)->'Point':
        if n==1:
            return point
        if n in memo:
            return memo[n]
        else:
            memo[n]=Point([newton(n-1).coordinates[0]-func.Map(newton(n-1))/diff(func,newton(n-1),0)])
            return memo[n]
    return newton(number)


    
#def integral(다중적분)
def integral(func:Function,bounds:List[Domain])->float:
    if func.dimention()!=len(bounds):
        raise ValueError("Function dimention does not match domain dimention")
    
    def integratevalue(d:int,nowpoint:Point)->float:
        if d==len(bounds):
            return func.Map(nowpoint)
        
        lower_bound=bounds[d].lowerbound.coordinates[0]
        upper_bound=bounds[d].upperbound.coordinates[0]

        result=0
        nowcoord=lower_bound
        while nowcoord<=upper_bound:
            nextpoint=Point(nowpoint.coordinates[:])
            nextpoint.coordinates[d]=nowcoord
            
            result+=integratevalue(d+1,nextpoint)*epsilon
            
            nowcoord+=epsilon
        return result
    
    initialpoint=Point([0]*func.dimention())
    return integratevalue(0,initialpoint)



#def line_integral(선적분)
def line_integral(func:Function,line:Parameter)->float:
    if not(func.dimention()==line.dimention()+1):
        raise ValueError("Dimention is not match")
    summ=0
    t=line.lowerbound.coordinates[0]
    tend=line.upperbound.coordinates[0]
    while t<=tend:
        nowpoint=line.Map(Point([t]))
        nextpoint=line.Map(Point([t+epsilon]))

        lenth=[
            (nextpoint.coordinates[i]-nowpoint.coordinates[i])/epsilon
            for i in range(line.dimention())
        ]


        ds=math.sqrt(sum(d**2 for d in lenth))
        funcvalue=func.Map(nowpoint)
        summ+=ds*funcvalue*epsilon
        t+=epsilon


    return summ


#def Graph2D(정의역 내 y=f(x) 그래프 그리기)
def Graph2D(func:Function,lowerbound:Point,upperbound:Point):

    xvalue=np.linspace(lowerbound.coordinates[0],upperbound.coordinates[0],1000)
    yvalues=[]

    for x in xvalue:
        try:
            y=func.Map(Point([x]))
            yvalues.append(y)
        except ValueError:
            yvalues.append(np.nan)

    plt.figure(figsize=(10,6))
    plt.plot(xvalue,yvalues,label="y=f(x)",color="blue")
    plt.axhline(0,color='black',linewidth=0.8,linestyle='--')
    plt.axvline(0,color='black',linewidth=0.8,linestyle='--')
    plt.title("Graph of y = f(x)")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.legend()
    plt.savefig("graph2D.png", dpi=300)
    plt.show()


#def solvediff(간단한 미분방정식 수치해석적 풀이)
def solvediff(diffequal: Callable[[float, float], float], y0: float, t_bounds: Point):
   
    t_start, t_end = t_bounds.coordinates
    t_values = [t_start]
    y_values = [y0]

    t = t_start
    y = y0

    while t < t_end:
        y += epsilon * diffequal(t, y)  
        t += epsilon
        t_values.append(t)
        y_values.append(y)

    plt.figure(figsize=(10, 6))
    plt.plot(t_values, y_values, label="y(t)", color="blue")
    plt.axhline(0, color='black', linewidth=0.8, linestyle='--')
    plt.axvline(0, color='black', linewidth=0.8, linestyle='--')
    plt.title("Solution of Differential Equation")
    plt.xlabel("t")
    plt.ylabel("y")
    plt.grid(True)
    plt.legend()
    plt.savefig("differential_solution_graph.png", dpi=300)
    plt.show()

    return [Point([t_values[i], y_values[i]]) for i in range(len(t_values))]


#testcase : class Point
point1=Point([1,1])
point2=Point([2,3])
k=2.3
print(point1.add(point2),point1.scale(k),point1.inner(point2),point1.degree(point2))


#testcase : class Function
def fx(point:Point)->float:
    x=point.coordinates[0]
    return x**2

def condition(x1:Point,x2:Point,point:Point)->bool:
    if x1.coordinates[0]<point.coordinates[0]<x2.coordinates[0]:
        return True
    else:
        return False

f_lb=Point([0])
f_ub=Point([5])
f_point=Point([2])
f=Function(f_lb,f_ub,fx,condition)
print(f.Map(f_point))


#testcase : class Parameter
p_lb=Point([0])
p_ub=Point([3])
p_point=Point([math.pi/2])
def parameter(point:Point)->'Point':
    t=point.coordinates[0]
    return Point([math.cos(t),math.sin(t)])
para=Parameter(p_lb,p_ub,parameter)
print(para.Map(p_point))


#testcase : def diff, def gradient
d_lb=Point([0,0])
d_ub=Point([10,10])
d_point=Point([2,4])
def d_condition(x1:Point,x2:Point,point:Point):
    for i in range(2):
        if not (x1.coordinates[i]<point.coordinates[i]<x2.coordinates[i]):
            return False
        else:
            return True
def fxy(point:Point)->float:
    x,y=point.coordinates
    return x**2+y**2
d_func=Function(d_lb,d_ub,fxy,d_condition)
print(diff(d_func,d_point,1),gradient(d_func,d_point))


#testcase : def newton-labson
n_lb=Point([0])
n_ub=Point([10])
n_point=Point([3])
def n_fx(point:Point)->float:
    x=point.coordinates[0]
    return x**2-2
n_func=Function(n_lb,n_ub,n_fx,condition)
print(newton_labson(n_func,n_point,5))


#testcase : def integral
i_lb=Point([0,0])
i_ub=Point([1,1])
def i_fxy(point:Point)->float:
    x,y=point.coordinates
    return x*y+x+y
def i_condition(x1:Point,x2:Point,point:Point)->bool:
    for i in range(2):
        if not (x1.coordinates[i]<=point.coordinates[i]<=x2.coordinates[i]):
            return False
    return True
i_func=Function(i_lb,i_ub,i_fxy,i_condition)
print(integral(i_func,[Domain(Point([0]),Point([1]),i_condition),Domain(Point([0]),Point([1]),i_condition)]))


#testcase : def line_integral
li_lb=Point([0,0])
li_ub=Point([1,1])
def li_fxy(point:Point)->float:
    x,y=point.coordinates
    return x**2+y**2
li_func=Function(li_lb,li_ub,li_fxy,i_condition)

l_lb=Point([0])
l_ub=Point([1])
def li_t(point:Point)->'Point':
    t=point.coordinates[0]
    return Point([t,t])
line=Parameter(l_lb,l_ub,li_t)

print(line_integral(li_func,line))



# testcase : def solvediff
def example_diffeq(t, y):
    return y  
y0 = 1
t_bounds = Point([-1, 6])  
solution = solvediff(example_diffeq, y0, t_bounds)
