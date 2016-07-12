import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
#from Bio import Entrez
from Bio import SeqIO
from Bio import Seq
from Bio.SeqUtils import GC

class MyWindow(Gtk.Window):
    
    def __init__(self):
        Gtk.Window.__init__(self)
        self.set_size_request(300, 100)
        self.button1_connection = 0
        
        self.set_title("Bio Programme")
        self.msg_label = Gtk.Label()
        self.board_label = Gtk.Label("Bio")
        self.entry = Gtk.Entry()
        self.button1 = Gtk.Button("GC content")
        self.button1.connect("clicked", self.clicked_callback)
        
        vbox = Gtk.VBox()
        bbox = Gtk.HButtonBox()
        bbox.add(self.button1)
        vbox.add(self.board_label)
        vbox.add(self.msg_label)
        vbox.add(self.entry)
        self.entry.connect("activate", self.enter_callback)
        vbox.add(bbox)
        self.add(vbox)
        
        self.show_all()
        self.connect("delete-event", self.on_quit)

    def on_quit(self, widget, event):
        Gtk.main_quit()

    def clicked_callback(self, button): #runs function when button is clicked
        self.gc_content()

    def enter_callback(self, widget): #runs function when enter pressed
        self.gc_content()
        
    def run(self):
        self.msg_label.set_text("What would you like to do?")
        self.set_button_labels("GC content")   
        self.enter_callback(Gtk.Entry, self.entry)
                
    def gc_content(self): #calculates gc % of sequence
        entry_text = self.entry.get_text()
        self.board_label.set_text("Type in your sequence")
        result = GC(entry_text)
        self.msg_label.set_text(str(result))


        

       
window = MyWindow()

Gtk.main()
