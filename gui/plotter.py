#! /usr/local/bin/python
from Tkinter import *
import math
import string
import os
import firstswig

root = Tk()

#first, a row for function entry and action button
fram = Frame(root)
Label(fram,text='dimensions').pack(side=LEFT)
ndims = Entry(fram)
ndims.insert(END, "2")
ndims.pack(side=LEFT, fill=BOTH, expand=1)
dbutt = Button(fram, text='Set')
dbutt.pack(side=LEFT)
Label(fram,text='f(x)').pack(side=LEFT)
func = Entry(fram)
func.pack(side=LEFT, fill=BOTH, expand=1)
butt = Button(fram, text='Plot')
butt.pack(side=RIGHT)
fram.pack(side=TOP)

#then a row to enter bound in
fram = Frame(root)
bounds = []
for label in 'minX','maxX','minY','maxY':
    Label(fram,text=label+':').pack(side=LEFT)
    edit = Entry(fram, width = 6)
    edit.pack(side=LEFT)
    bounds.append(edit)
fram.pack(side=TOP)
fram.pack_forget()
fram.pack(side=TOP)

dimframes = []
mlabels = []
minvals = []
maxvals = []
radiobuttons = []
pyPointDef = firstswig.doubleArray(1)
pyPointDefLength = 1

x = 2

lframe = Frame(root)
lframe.pack(side = TOP)
v = IntVar()

def fillValues():
  global v, maxvals, minvals, pyPointDef, pyPointDefLength
  for i in range(len(minvals)):
    minvals[i].delete(0,END)
    maxvals[i].delete(0,END)
    if i * 2 + 1 >= pyPointDefLength:
      minvals[i].insert(END, pyPointDef[i*2])
      if v.get() == i: 
        maxvals[i].insert(END, pyPointDef[i*2+1])
    else:
      minvals[i].insert(END, -1.0)
      if v.get() == i: 
        maxvals[i].insert(END, 1.0)

  
def radioSelection():
  global v, mlabels, maxvals
  for i in range(x):
    if v.get() == i:
      mlabels[i].config(text="Range:")
      maxvals[i].config(borderwidth=2, highlightthickness=1,
        takefocus=True, state=NORMAL)
      if i*2+1 < pyPointDefLength:
        maxvals[i].delete(0,END)
        maxvals[i].insert(END,pyPointDef[i*2+1])
    else:
      mlabels[i].config(text="Value:")
      maxvals[i].delete(0,END)
      maxvals[i].config(borderwidth=0,highlightthickness=0,
        takefocus=False, state=DISABLED)

def readValues():
  global pyPointDef, pyPointDefLength, minvals, maxvals
  pyPointDefLength = len(minvals)*2
  pyPointDef = firstswig.doubleArray(pyPointDefLength)
  for i in range(len(minvals)):
    try: pyPointDef[i*2] = float(minvals[i].get()) 
    except: pyPointDef[i*2] = 0.0
    try: pyPointDef[i*2+1] = float(maxvals[i].get()) 
    except: pyPointDef[i*2+1] = 0.0

def removeFrames():
  global dimframes, mlabels, minvals, maxvals, radiobuttons
  readValues()
  for ml in mlabels: ml.destroy()
  for mi in minvals: mi.destroy()
  for ma in maxvals: ma.destroy()
  for ra in radiobuttons: ra.destroy() 
  for l in dimframes: l.destroy()
  dimframes = []
  mlabels = []
  minvals = []
  maxvals = []
  radiobuttons = []
  

def addFrames():
    global x, mlabels
    removeFrames()
    x = int(ndims.get())
    if x < 0: x = 0
    for i in range(x):
      dimframes.append( Frame(lframe) )
      mlabels.append( Label(dimframes[i], text = 'Value: '))
      mlabels[i].pack(side=LEFT)
      minvals.append( Entry(dimframes[i], width=5) )
      minvals[i].pack(side=LEFT)
      maxvals.append( Entry(dimframes[i], width=5))
      maxvals[i].pack(side=LEFT)
      radiobuttons.append( Radiobutton(dimframes[i],variable=v,value=i,
	 command=radioSelection))
      radiobuttons[i].pack()
      dimframes[i].pack()
      v.set(0)
      fillValues()
    radioSelection()
addFrames()
radioSelection()
    
    

#and finally the canvas
c = Canvas(root)
c.pack(side=TOP, fill=BOTH, expand=1)

cedit = Entry(root)
cedit.pack(side=TOP, fill=BOTH, expand=1)

myobj = firstswig.FirstClass(4)
k = 1000
def enterPressedCommand(*ignore):
    global myobj
    global k
    str = cedit.get()
    print str
    t.insert(END, str + "\n")
    cedit.delete(0,END)
    sed_command = "cat template | sed 's/\/\/REPLACE_ME/" + str \
        + "/' > firstswig.cxx"
    print sed_command
    os.system(sed_command)
    os.system("./mymake")
    #import firstswig
    #reload(firstswig._firstswig)
    reload(firstswig)
    
    myobj = firstswig.FirstClass(k)
    k = k*2
    myobj.printVal("Will it work?")
    
    myarray =  firstswig.doubleArray(3)
    myarray[0] = 1.0
    myarray[1] = 2.0
    myarray[2] = 4.0
    
    #x = myobj.shiftArray(myarray, 3)
    #print x
    #x = myobj.shiftArray(myarray, 3)
    #print x
    y= myobj.myevaluate(myarray, 3)
    print y

    

cedit.bind('<Return>',enterPressedCommand)

t = Text(root)
t.pack(side=TOP, fill=BOTH, expand=1)
t.start_command = '1.0'
def enterPressed(*ignore): 
    s = t.get(t.start_command, END)
    print "[", string.strip(s), "]"
    print t.start_command
    t.insert(END, "\n_____ ")
    t.start_command = CURRENT

t.bind('<Return>', enterPressed)

def minimax(values=[0.0, 1.0, 0.0, 1.0]):
    "Adjust and display X and Y bounds"
    for i in range(4):
        edit = bounds[i]
        try: values[i] = float(edit.get())
        except: pass
        edit.delete(0, END)
        edit.insert(END, '%.2f'%values[i])
    return values

def plot():
    "Plot given functions with given bounds"
    minx, maxx, miny, maxy = minimax()
    
    #get and compile the function
    f = func.get()
    f = compile(f, f, 'eval')
    
    #get Canvas X and Y dimensions
    CX = c.winfo_width()
    CY = c.winfo_height()

    #compute coordinates for line
    coords = []
    for i in range(0,CX,5):
        coords.append(i)
        x = minx + ((maxx-minx)*i)/CX
        z = 5
        y = eval(f, vars(math), {'x':x,'z':z})
        j = CY - CY*(y-miny)/(maxy - miny)
        coords.append(j)

    # draw line
    c.delete(ALL)
    c.create_line(*coords)
 
butt.config(command=plot)
dbutt.config(command=addFrames)
minvals[0].config(state=DISABLED)


# give an initial example in lieu of docs
f = 'sin(x) + cos(x)'
func.insert(END, f)
minimax([0.0, 10.0, -2.0, 2.0])

cedit.delete(0,END)
cedit.insert(END, "result = vals[0] + vals[1];")
cedit.focus_set()
root.mainloop()
