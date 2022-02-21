
from tkinter import *
from tkinter.ttk import Frame, Button, Scrollbar

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import numpy as np

from abc import ABC, abstractmethod


class Model(ABC):
    def __init__(self):
        self._observers = list()
        self.views = self._get_views()
        for view in self.views:
            self.attach(view)

    def attach(self, observer):
        self._observers.append(observer)

    def notify(self):
        for observer in self._observers:
            observer.update()

    @abstractmethod
    def _get_views(self):
        raise NotImplemented()

    def controls(self, container):
        return list()

class View(ABC):
    def __init__(self, model):
        self.model = model
        self._active = False
        self.viewport = None

    def show(self, viewport):
        self.viewport = viewport
        self.update()

    def update(self):
        if not self.viewport is None:
            figure = self.refresh()
            for widget in self.viewport.winfo_children():
                widget.destroy()

            chart_type = FigureCanvasTkAgg(figure, self.viewport)
            chart_type.get_tk_widget().pack(fill = BOTH, expand = True)

    @abstractmethod
    def refresh(self):
        raise NotImplemented()

    def view_controls(self):
        return list()
   
class Slider(Frame):

    def __init__(self, container, name, low, high, value, command = None):
        super().__init__(container)

        self.title_frame = Frame(self)
        self.title_frame.pack(side = 'top', fill = 'x')
        # self.slider_frame = Frame(self)
        # self.slider_frame.pack(side = 'top', fill = 'x')


        name_label = Label(self.title_frame,
                           text = name,
                           width = 10,
                           justify = LEFT)
        name_label.pack(side = 'left')

        from_value = low
        to_value = high

        low_label = Label(self.title_frame,
                          text = f'{low:.2g}')
        low_label.pack(side = 'left')
        self.scale = Scale(self.title_frame,
                          orient=HORIZONTAL,
                          length=200,
                          width=10,
                          resolution = (to_value - from_value)/1000,
                          sliderlength=10,
                          from_=from_value,
                          to=to_value,
                          command = command)
        self.scale.pack(side = 'left', fill = 'x')
        self.scale.set(value)

        high_label = Label(self.title_frame,
                           text = f'{high:.2g}')
        high_label.pack(side = 'left')


    def get(self):
        return self.scale.get()


class Window():
    def __init__(self, master, model):
        self.master = master
        self.master.geometry("1200x800+300+300")

        self.model = model

        self.view_controls()

        self.control_frame = Frame(relief=RAISED, borderwidth=1)
        self.control_frame.pack(side = 'left', fill = 'y')
        
        #scrollbar = Scrollbar(self.control_frame, orient = 'vertical')
        #scrollbar.pack(side = RIGHT, fill = 'y')

        controls = self.model.controls(self.control_frame)

        if isinstance(controls, dict):
            controls = list(controls.values())
        for ctr in controls:
            ctr.pack()

        self.view_frame = Frame(relief=RAISED, borderwidth=1)
        self.view_frame.pack(side = 'left', fill = BOTH, expand = True)

        self._active_view_index = -1
        self.next_view()

    def view_controls(self):
        self.view_frame = Frame(relief=RAISED, borderwidth=1)
        self.view_frame.pack(side = 'top', fill = 'x')

        prev = Button(self.view_frame,
                      text = 'prev',
                      command = self.prev_view,
                      width = 20)
        prev.pack(side = 'left')
        next = Button(self.view_frame,
                      text = 'next',
                      command = self.next_view,
                      width = 20)
        next.pack(side = 'left')


    def next_view(self):
        if len(self.model.views) > 0:
            self.model.views[self._active_view_index].viewport = None
            self._active_view_index = ((self._active_view_index+1) % len(self.model.views))
            self.model.views[self._active_view_index].show(self.view_frame)

    def prev_view(self):
        if len(self.model.views) > 0:
            self.model.views[self._active_view_index].viewport = None
            self._active_view_index = ((self._active_view_index-1) % len(self.model.views))
            self.model.views[self._active_view_index].show(self.view_frame)


class MockView1(View):
    def refresh(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot([0,1],self.model.data, self.model.col)
        plt.title('test plot')
        plt.xlim([-1,2])
        plt.ylim([0,1])
        return fig

class MockView2(View):
    def refresh(self):
        fig = plt.figure()
        plt.plot([0,1], self.model.data, 'go-')
        plt.title('test plot')
        plt.xlim([-1,2])
        plt.ylim([0,1])
        return fig


class MockModel(Model):

    def __init__(self):
        self.col = 'r-'
        self.data = np.array([.1,.8])
        super().__init__()

    def _get_views(self):
        view1 = MockView1(self)
        view2 = MockView2(self)
        return [view1, view2]

    def change_color(self):
        self.col = 'b-' if self.col == 'r-' else 'r-'
        self.notify()

    def random_data(self):
        self.data = np.random.uniform(0,1,(2,))
        self.notify()

    def controls(self, container):
        contr = [
                Button(container,
                       text = 'col',
                       command = self.change_color),
                Button(container,
                       text = 'data',
                       command = self.random_data)
                 ]


        return contr
# choice:
#  pick geometry
#   change geometry
#  pick segment (if geometry exists)
#    change boundary type
#  define global body force
#  define global material
#
# switch view: solution (components), error, predicted h
# side-by-side view of predicted and adaptive
#
# show individual contributions (kernels, corners, responsibilities)
#
# field to show stats (time taken, error, ...)
#

def run_demo(model):

    # GUI init
    root = Tk()
    #root.tk.call('tk', 'scaling', 0.7)
    root.bind('<Escape>', lambda _: root.quit())
    d = Window(root, model)
    root.mainloop()

    # pick polygon from canvas
    # using sliders:
    # specify bc
    # specify material
    # specify body force
    #

    # grab geometry/physics, build bvp

    # predict mesh

    # solve
    #


if __name__ == '__main__':
    model = MockModel()
    run_demo(model)
