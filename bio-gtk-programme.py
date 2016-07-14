import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio.Alphabet import IUPAC


class LnFrame(Gtk.Bin):
    def __init__(self):
        Gtk.Bin.__init__(self)

        self.msg_label = Gtk.Label()
        self.board_label = Gtk.Label("Enter sequence")
        self.entry = Gtk.Entry()
        self.button = Gtk.Button("Sequence Length")
        self.button.connect("clicked", self.clicked_callback)
        bbox = Gtk.HButtonBox()
        bbox.add(self.button)
        ln_box = Gtk.VBox()
        ln_box.add(self.board_label)
        ln_box.add(self.msg_label)
        ln_box.add(self.entry)
        self.entry.connect("activate", self.enter_callback)
        ln_box.add(bbox)
        self.add(ln_box)

    def seq_len(self): #gives the length of the sequence
        entry_text = Seq(self.entry.get_text(), IUPAC.unambiguous_dna)
        len_result = len(entry_text)
        self.board_label.set_text("Result")
        self.msg_label.set_text(str(len_result))

    def clicked_callback(self, button): #runs function when button is clicked
        self.seq_len()
        
    def enter_callback(self, widget): #runs function when enter pressed
        self.seq_len()
        

class CmFrame(Gtk.Bin):
    def __init__(self):
        Gtk.Bin.__init__(self)

        self.msg_label = Gtk.Label()
        self.board_label = Gtk.Label("Enter sequences")
        self.entry1 = Gtk.Entry()
        self.entry2 = Gtk.Entry()
        self.button = Gtk.Button("Sequence Comparison")
        self.button.connect("clicked", self.clicked_callback)
        bbox = Gtk.HButtonBox()
        bbox.add(self.button)
        cm_box = Gtk.VBox()
        cm_box.add(self.board_label)
        cm_box.add(self.msg_label)
        cm_box.add(self.entry1)
        cm_box.add(self.entry2)
        self.entry1.connect("activate", self.enter_callback)
        self.entry2.connect("activate", self.enter_callback)
        cm_box.add(bbox)
        self.add(cm_box)

    def seq_comp(self): #compares two sequences
        entry_text1 = Seq(self.entry1.get_text(), IUPAC.unambiguous_dna)
        entry_text2 = Seq(self.entry2.get_text(), IUPAC.unambiguous_dna)
        comp_result = ((entry_text1) == (entry_text2))
        self.board_label.set_text("Result")
        self.msg_label.set_text(str(comp_result))

    def clicked_callback(self, button): #runs function when button is clicked
        self.seq_comp()
        
    def enter_callback(self, widget): #runs function when enter pressed
        self.seq_comp()
        

class GcFrame(Gtk.Bin):
    def __init__(self):
        Gtk.Bin.__init__(self)

        self.msg_label = Gtk.Label()
        self.board_label = Gtk.Label("Enter sequence")
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
        gc_result = GC(entry_text)
        self.board_label.set_text("Result")
        self.msg_label.set_text(str(gc_result))

    def clicked_callback(self, button): #runs function when button is clicked
        self.gc_content()
        
    def enter_callback(self, widget): #runs function when enter pressed
        self.gc_content()
        
class RcFrame(Gtk.Bin):
    def __init__(self):
        Gtk.Bin.__init__(self)

        self.msg_label = Gtk.Label()
        self.board_label = Gtk.Label("Enter sequence")
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
        rc_result = entry_text.reverse_complement()
        self.board_label.set_text("Result")
        self.msg_label.set_text(str(rc_result))

    def clicked_callback(self, button): #runs function when button is clicked
        self.rev_comp()
        
    def enter_callback(self, widget): #runs function when enter pressed
        self.rev_comp()


class TrFrame(Gtk.Bin):
    def __init__(self):
        Gtk.Bin.__init__(self)

        self.msg_label = Gtk.Label()
        self.board_label = Gtk.Label("Enter sequence")
        self.entry = Gtk.Entry()
        self.button = Gtk.Button("Translate Sequence")
        self.button.connect("clicked", self.clicked_callback)
        bbox = Gtk.HButtonBox()
        bbox.add(self.button)
        tr_box = Gtk.VBox()
        tr_box.add(self.board_label)
        tr_box.add(self.msg_label)
        tr_box.add(self.entry)
        self.entry.connect("activate", self.enter_callback)
        tr_box.add(bbox)
        self.add(tr_box)

    def translation(self): #translates sequence into protein
        entry_text = Seq(self.entry.get_text(), IUPAC.unambiguous_dna)
        mrna = entry_text.transcribe()
        tr_result = mrna.translate()
        self.board_label.set_text("Result")
        self.msg_label.set_text(str(tr_result))

    def clicked_callback(self, button): #runs function when button is clicked
        self.translation()
        
    def enter_callback(self, widget): #runs function when enter pressed
        self.translation()


class EnFrame(Gtk.Bin):
    def __init__(self):
        Gtk.Bin.__init__(self)

        self.msg_label = Gtk.Label()
        self.board_label = Gtk.Label("Enter ID")
        self.entry = Gtk.Entry()
        self.button = Gtk.Button("Entrez")
        self.button.connect("clicked", self.clicked_callback)
        bbox = Gtk.HButtonBox()
        #fbox = Gtk.FileChooserButton()
        #fbox.connect("file_set", self.on_file_selected)        
        bbox.add(self.button)
        en_box = Gtk.VBox()
        en_box.add(self.board_label)
        en_box.add(self.msg_label)
        en_box.add(self.entry)
        self.entry.connect("activate", self.enter_callback)
        en_box.add(bbox)
        #en_box.add(fbox)
        self.add(en_box)

    def on_file_selected(self, entry):
        print("You chose ", entry.get_filename())
        
    def entrez_db(self): #finds details from Entrez database
        entry_text = Entrez.efetch(db = "nucleotide", id = [self.entry.get_text()], rettype = "fasta")
        en_result = list(SeqIO.parse(entry_text, "fasta"))
        self.board_label.set_text("Result")
        self.msg_label.set_text(str("Written to file 'bio-gtk-entrez.txt'"))#en_result))
        SeqIO.write(en_result, "bio-gtk-entrez.txt", "fasta")

    def clicked_callback(self, button): #runs function when button is clicked
        self.entrez_db()
        
    def enter_callback(self, widget): #runs function when enter pressed
        self.entrez_db()


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

        ln_box = LnFrame()
        cm_box = CmFrame()
        gc_box = GcFrame()
        rc_box = RcFrame()
        tr_box = TrFrame()
        en_box = EnFrame()

        stack.add_titled(ln_box, "ln", "Sequence Length")
        stack.add_titled(cm_box, "cm", "Sequence Comparison")
        stack.add_titled(gc_box, "gc", "Calculate GC")
        stack.add_titled(rc_box, "rc", "Reverse Complement")
        stack.add_titled(tr_box, "tr", "Translate Sequence")
        stack.add_titled(en_box, "en", "Entrez")
        self.add(vbox)
        
        self.show_all()
        self.connect("delete-event", self.on_quit)
        

    def on_quit(self, widget, event):
        Gtk.main_quit()

       
window = MyWindow()

Gtk.main()
