import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
#from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio.Alphabet import IUPAC


class GcFrame(Gtk.Bin):
    def __init__(self):
        Gtk.Bin.__init__(self)

        self.msg_label = Gtk.Label()
        self.board_label = Gtk.Label("Bio")
        self.entry = Gtk.Entry()
        self.button = Gtk.Button("GC content")
        self.button.connect("clicked", self.clicked_callback)
        bbox = Gtk.HButtonBox()
        bbox.add(self.button)
        gc_box = Gtk.VBox()
        gc_box.add(self.board_label)
        gc_box.add(self.msg_label)
        gc_box.add(self.entry)
        self.entry.connect("activate", self.enter_callback)
        gc_box.add(bbox)
        self.add(gc_box)

    def gc_content(self): #calculates gc % of sequence
        entry_text = Seq(self.entry.get_text(), IUPAC.unambiguous_dna)
        self.board_label.set_text("Type in your sequence")
        gc_result = GC(entry_text)
        self.msg_label.set_text(str(gc_result))

    def clicked_callback(self, button): #runs function when button is clicked
        self.gc_content()
        
    def enter_callback(self, widget): #runs function when enter pressed
        self.gc_content()
        
class RcFrame(Gtk.Bin):
    def __init__(self):
        Gtk.Bin.__init__(self)

        self.msg_label = Gtk.Label()
        self.board_label = Gtk.Label("Bio")
        self.entry = Gtk.Entry()
        self.button = Gtk.Button("Reverse Complement")
        self.button.connect("clicked", self.clicked_callback)
        bbox = Gtk.HButtonBox()
        bbox.add(self.button)
        rc_box = Gtk.VBox()
        rc_box.add(self.board_label)
        rc_box.add(self.msg_label)
        rc_box.add(self.entry)
        self.entry.connect("activate", self.enter_callback)
        rc_box.add(bbox)
        self.add(rc_box)

    def rev_comp(self): #gives the reverse complement of DNA sequence
        entry_text = Seq(self.entry.get_text(), IUPAC.unambiguous_dna)
        self.board_label.set_text("Type in your sequence")
        rc_result = entry_text.reverse_complement()
        self.msg_label.set_text(str(rc_result))

    def clicked_callback(self, button): #runs function when button is clicked
        self.rev_comp()
        
    def enter_callback(self, widget): #runs function when enter pressed
        self.rev_comp()


class MyWindow(Gtk.Window):
    
    def __init__(self):
        Gtk.Window.__init__(self)
        self.set_size_request(300, 100)
        
        self.set_title("Bio Programme")
        self.msg_label = Gtk.Label()
        self.board_label = Gtk.Label("Bio")
        self.entry = Gtk.Entry()
        
        vbox = Gtk.VBox()
        stack = Gtk.Stack()
        stack.set_transition_type(Gtk.StackTransitionType.SLIDE_LEFT_RIGHT)
        stack.set_transition_duration(250)
        switcher = Gtk.StackSwitcher()
        switcher.set_stack(stack)
        vbox.add(switcher)
        vbox.add(stack)

        gc_box = GcFrame()
        rc_box = RcFrame()
        
        stack.add_titled(gc_box, "gc", "Calculate GC")
        stack.add_titled(rc_box, "rc", "Reverse Complement")
        
        self.add(vbox)
        
        self.show_all()
        self.connect("delete-event", self.on_quit)
        

    def on_quit(self, widget, event):
        Gtk.main_quit()

       
window = MyWindow()

Gtk.main()
