
# coding: utf-8

# In[1]:


range(5, 10)


# In[2]:


2+2


# In[3]:


for n in range(2, 10):
    for x in range(2, n):
        if n % x == 0:
            print(n, 'equals', x, '*', n//x)
            break
    else:
        # loop fell through without finding a factor
        print(n, 'is a prime number')


# In[4]:


def fib(n):    # write Fibonacci series up to n
    """Print a Fibonacci series up to n."""
    a, b = 0, 1
    while a < n:
        print(a, end=' ')
        a, b = b, a+b
    print()

# Now call the function we just defined:
fib(2000)


# In[5]:


def fib2(n):  # return Fibonacci series up to n
    """Return a list containing the Fibonacci series up to n."""
    result = []
    a, b = 0, 1
    while a < n:
        result.append(a)    # see below
        a, b = b, a+b
    return result

f100 = fib2(100)    # call it
f100                # write the result


# In[6]:


squares = []
for x in range(10):
    squares.append(x**2)

squares


# In[7]:


matrix = [
    [1, 2, 3, 4],
    [5, 6, 7, 8],
    [9, 10, 11, 12],
]


# In[8]:


[[row[i] for row in matrix] for i in range(4)]


# In[9]:


transposed = []
for i in range(4):
    # the following 3 lines implement the nested listcomp
    transposed_row = []
    for row in matrix:
        transposed_row.append(row[i])
    transposed.append(transposed_row)

transposed


# In[10]:


a = [-1, 1, 66.25, 333, 333, 1234.5]
del a[0]
a
[1, 66.25, 333, 333, 1234.5]
del a[2:4]
a
[1, 66.25, 1234.5]
del a[:]
a


# In[12]:


# Fibonacci numbers module

def fib(n):    # write Fibonacci series up to n
    a, b = 0, 1
    while b < n:
        print(b, end=' ')
        a, b = b, a+b
    print()

def fib2(n):   # return Fibonacci series up to n
    result = []
    a, b = 0, 1
    while b < n:
        result.append(b)
        a, b = b, a+b
    return result


# In[13]:


import fibo


# In[14]:


import math
print('The value of PI is approximately {0:.3f}.'.format(math.pi))


# In[15]:


with open('http://isbweb.org/data/pezzack/Pezzack.txt') as f:
    read_data = f.read()
f.closed


# In[16]:


open('http://isbweb.org/data/pezzack/Pezzack.txt','r+')


# In[17]:


f = open('http://isbweb.org/data/pezzack/Pezzack.txt', 'w')


# In[19]:


f=open(Users\caiol\Desktop\exerc_phyton.txt)


# In[20]:


f=open(http://localhost:8888/notebooks/Desktop/exerc_phyton.txt)

