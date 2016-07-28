#CRASHES WHEN TRYING TO GET LARGE FILES FROM ENTREZ
#can't deal with N in translation
#don't know how to access data from list store

import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
from Bio import Entrez, SeqIO
#Entrez.email = "A.N.Other@example.com" 
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio.Alphabet import IUPAC


class LnFrame(Gtk.Bin):
    def __init__(self, open_sequences):
        Gtk.Bin.__init__(self)

        self.open_sequences = open_sequences
        
        self.builder = Gtk.Builder()
        self.builder.add_from_file("seqlen-page-glade.glade")
        self.sl_box = self.builder.get_object("SeqLen-box")
        self.add(self.sl_box)
    
        #self.entry = self.builder.get_object("SeqLen-entry")
        self.cbox = self.builder.get_object("SeqLen-cbox")
        self.c_entry = self.builder.get_object("SeqLen-cbox-entry")
        self.button = self.builder.get_object("SeqLen-button")
        self.label = self.builder.get_object("SeqLen-result")
        self.button.connect("clicked", self.clicked_callback)
        self.c_entry.connect("activate", self.enter_callback)

        renderer = Gtk.CellRendererText()
        self.cbox.pack_start(renderer, True)
        self.cbox.add_attribute(renderer, "text", 0)
        self.cbox.set_model(open_sequences)
        
    def seq_len(self): #gives the length of the sequence
        entry_text = Seq(self.c_entry.get_text(), IUPAC.unambiguous_dna)
        len_result = len(entry_text)
        self.label.set_text(str(len_result))

    def clicked_callback(self, button): #runs function when button is clicked
        self.seq_len()
        
    def enter_callback(self, widget): #runs function when enter pressed
        self.seq_len()
        #self.open_sequences.append([self.c_entry.get_text()]) #adds typed text to liststore

        

class CmFrame(Gtk.Bin):
    def __init__(self, open_sequences):
        Gtk.Bin.__init__(self)

        self.open_sequences = open_sequences

        self.builder = Gtk.Builder()
        self.builder.add_from_file("seqcomp-page-glade.glade")
        self.cm_box = self.builder.get_object("SeqComp-box")
        self.add(self.cm_box)

        #self.entry1 = self.builder.get_object("SeqComp-entry1")
        #self.entry2 = self.builder.get_object("SeqComp-entry2")
        self.cbox1 = self.builder.get_object("SeqComp-cbox1")
        self.c_entry1 = self.builder.get_object("SeqComp-cbox-entry1")
        self.cbox2 = self.builder.get_object("SeqComp-cbox2")
        self.c_entry2 = self.builder.get_object("SeqComp-cbox-entry2")
        self.button = self.builder.get_object("SeqComp-button")
        self.label = self.builder.get_object("SeqComp-result")
        self.button.connect("clicked", self.clicked_callback)
        
        self.c_entry1.connect("activate", self.enter_callback)
        self.c_entry2.connect("activate", self.enter_callback)

        renderer = Gtk.CellRendererText()
        self.cbox1.pack_start(renderer, True)
        self.cbox1.add_attribute(renderer, "text", 0)
        self.cbox1.set_model(open_sequences)
        self.cbox2.pack_start(renderer, True)
        self.cbox2.add_attribute(renderer, "text", 0)
        self.cbox2.set_model(open_sequences)

    def seq_comp(self): #compares two sequences
        entry_text1 = Seq(self.c_entry1.get_text(), IUPAC.unambiguous_dna)
        entry_text2 = Seq(self.c_entry2.get_text(), IUPAC.unambiguous_dna)
        comp_result = ((entry_text1) == (entry_text2))
        self.label.set_text(str(comp_result))

    def clicked_callback(self, button): #runs function when button is clicked
        self.seq_comp()
        
    def enter_callback(self, widget): #runs function when enter pressed
        self.seq_comp()
        

class GcFrame(Gtk.Bin):
    def __init__(self, open_sequences):
        Gtk.Bin.__init__(self)

        self.open_sequences = open_sequences

        self.builder = Gtk.Builder()
        self.builder.add_from_file("gc-page-glade.glade")
        self.gc_box = self.builder.get_object("GC-box")
        self.add(self.gc_box)

        self.cbox = self.builder.get_object("GC-cbox")
        self.c_entry = self.builder.get_object("GC-cbox-entry")
        self.button = self.builder.get_object("GC-button")
        self.label = self.builder.get_object("GC-result")
        self.button.connect("clicked", self.clicked_callback)
        self.c_entry.connect("activate", self.enter_callback)
        
        renderer = Gtk.CellRendererText()
        self.cbox.pack_start(renderer, True)
        self.cbox.add_attribute(renderer, "text", 0)
        self.cbox.set_model(open_sequences)
        
    def gc_content(self): #calculates gc % of sequence
        entry_text = Seq(self.c_entry.get_text(), IUPAC.unambiguous_dna)
        gc_result = GC(entry_text)
        self.label.set_text(str(gc_result))

    def clicked_callback(self, button): #runs function when button is clicked
        self.gc_content()
        
    def enter_callback(self, widget): #runs function when enter pressed
        self.gc_content()
        

class RcFrame(Gtk.Bin):
    def __init__(self, open_sequences):
        Gtk.Bin.__init__(self)

        self.open_sequences = open_sequences
        
        self.builder = Gtk.Builder()
        self.builder.add_from_file("revcomp-page-glade.glade")
        self.rc_box = self.builder.get_object("RevComp-box")
        self.add(self.rc_box)

        self.cbox = self.builder.get_object("RevComp-cbox")
        self.c_entry = self.builder.get_object("RevComp-cbox-entry")
        self.button = self.builder.get_object("RevComp-button")
        self.label = self.builder.get_object("RevComp-result")
        self.button.connect("clicked", self.clicked_callback)
        self.c_entry.connect("activate", self.enter_callback)

        renderer = Gtk.CellRendererText()
        self.cbox.pack_start(renderer, True)
        self.cbox.add_attribute(renderer, "text", 0)
        self.cbox.set_model(open_sequences)

    
    def rev_comp(self): #gives the reverse complement of DNA sequence
        entry_text = Seq(self.c_entry.get_text())
        rc_result = entry_text.reverse_complement()
        if len(rc_result) < 30:
            self.label.set_text(str(rc_result))
        else:
            tv = Gtk.TextView()
            tv.get_buffer().set_text(str(rc_result))
            tv.set_editable(False)          
            sw = Gtk.ScrolledWindow()
            sw.set_size_request(300,200)
            sw.add(tv)
            w = Gtk.Window()
            w.add(sw)
            w.show_all()

    def clicked_callback(self, button): #runs function when button is clicked
        self.rev_comp()
        
    def enter_callback(self, widget): #runs function when enter pressed
        self.rev_comp()


class TrFrame(Gtk.Bin):
    def __init__(self, open_sequences):
        Gtk.Bin.__init__(self)

        self.open_sequences = open_sequences

        self.builder = Gtk.Builder()
        self.builder.add_from_file("transl-page-glade.glade")
        self.tr_box = self.builder.get_object("Transl-box")
        self.add(self.tr_box)

        self.cbox = self.builder.get_object("Transl-cbox")
        self.c_entry = self.builder.get_object("Transl-cbox-entry")
        self.button1 = self.builder.get_object("Transl-button")
        self.button2 = self.builder.get_object("Codon-button")
        self.label = self.builder.get_object("Transl-result")
        self.button1.connect("clicked", self.clicked_callback1)
        self.c_entry.connect("activate", self.enter_callback)
        self.button2.connect("clicked", self.clicked_callback2)

        renderer = Gtk.CellRendererText()
        self.cbox.pack_start(renderer, True)
        self.cbox.add_attribute(renderer, "text", 0)
        self.cbox.set_model(open_sequences)
        
    def translation(self): #translates sequence into protein
        entry_text = Seq(self.c_entry.get_text(), IUPAC.unambiguous_dna)
        mrna = entry_text.transcribe()
        tr_result = mrna.translate()
        if len(tr_result) < 100:
            self.label.set_text(str(tr_result))
        else:
            tv = Gtk.TextView()
            tv.get_buffer().set_text(str(tr_result))
            tv.set_editable(False)          
            sw = Gtk.ScrolledWindow()
            sw.set_size_request(300,200)
            sw.add(tv)
            w = Gtk.Window()
            w.add(sw)
            w.show_all()

    def clicked_callback1(self, button1): #runs function when button is clicked
        self.translation()

    def codon(self):
        from Bio.Data import CodonTable
        standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
        self.label.set_text(str(standard_table))

    def clicked_callback2(self, button2):
        self.codon()
                            
    def enter_callback(self, widget): #runs function when enter pressed
        self.translation()


class EnFrame(Gtk.Bin):
    def __init__(self, open_sequences):
        Gtk.Bin.__init__(self)

        self.open_sequences = open_sequences

        self.builder = Gtk.Builder()
        self.builder.add_from_file("entrez-page-glade.glade")
        self.en_box = self.builder.get_object("Entrez-box")
        self.add(self.en_box)

        self.cbox = self.builder.get_object("Entrez-cbox")
        self.c_entry = self.builder.get_object("Entrez-cbox-entry")
        self.button = self.builder.get_object("Entrez-button")
        self.label_e = self.builder.get_object("Entrez-label")
        self.label = self.builder.get_object("Entrez-result")
        self.label1 = self.builder.get_object("Entrez-result1")
        self.button.connect("clicked", self.clicked_callback)
        self.fbox = self.builder.get_object("Entrez-file")
        self.fbox.connect("file_set", self.on_file_selected)        
        self.c_entry.connect("activate", self.enter_callback)

        renderer = Gtk.CellRendererText()
        self.cbox.pack_start(renderer, True)
        self.cbox.add_attribute(renderer, "text", 0)
        self.cbox.set_model(open_sequences)

    def on_file_selected(self, entry): #opens file in window
        file = open(self.fbox.get_filename())
        tv = Gtk.TextView()
        tv.get_buffer().set_text(file.read())
        tv.set_editable(False)
        sw = Gtk.ScrolledWindow()
        sw.set_size_request(400, 400)
        sw.add(tv)
        w = Gtk.Window()
        w.add(sw)
        w.show_all()
        self.open_sequences.append([self.fbox.get_filename()]) #adds filename to liststore
    
        
    def entrez_db(self): #finds details from Entrez database
        entry_text = Entrez.efetch(db = "nucleotide", id = [self.c_entry.get_text()], rettype = "fasta")
        en_result = SeqIO.read(entry_text, "fasta")
        self.label_e.set_text(str("Written to file (Entrez ID)"))
        SeqIO.write(en_result, en_result.id, "fasta")
        self.label.set_text(str(en_result.description))
        self.label1.set_text(str(en_result.seq))       
        
    def clicked_callback(self, button): #runs function when button is clicked
        self.entrez_db()
        
    def enter_callback(self, widget): #runs function when enter pressed
        self.open_sequences.append([self.c_entry.get_text()]) #adds typed text to liststore



class OsFrame(Gtk.Bin): #opens sequences for later use
    def __init__(self, open_sequences):

        Gtk.Bin.__init__(self)

        self.open_sequences = open_sequences

        self.builder = Gtk.Builder()
        self.builder.add_from_file("openseq-page.glade")
        self.os_box = self.builder.get_object("Open-box")
        self.add(self.os_box)

        self.cbox = self.builder.get_object("Open-cbox")
        self.c_entry = self.builder.get_object("Open-cbox-entry")
        self.fbox = self.builder.get_object("Open-file")
        self.fbox.connect("file_set", self.on_file_selected)        
        self.c_entry.connect("activate", self.enter_callback)

        renderer = Gtk.CellRendererText()
        self.cbox.pack_start(renderer, True)
        self.cbox.add_attribute(renderer, "text", 0)
        self.cbox.set_model(open_sequences)
        
    def on_file_selected(self, entry): #opens a file and adds it to list store
        file = open(self.fbox.get_filename())
        self.open_sequences.append([self.fbox.get_filename()])

    def enter_callback(self, widget): #runs function when enter pressed
        self.open_sequences.append([self.c_entry.get_text()]) #adds typed text to liststore
        
        
class MyWindow(Gtk.Window):
    
    def __init__(self):
        Gtk.Window.__init__(self)
        self.set_size_request(300, 300)
        
        vbox = Gtk.VBox()       
        nbook = Gtk.Notebook()
        nbook.set_tab_pos(Gtk.PositionType.LEFT)
        nbook.set_scrollable(True)
        vbox.add(nbook)

        self.open_sequences = Gtk.ListStore(str)

        self.set_title("Bio Programme")
                          
        os_box = OsFrame(self.open_sequences)
        en_box = EnFrame(self.open_sequences)
        sl_box = LnFrame(self.open_sequences)
        cm_box = CmFrame(self.open_sequences)
        gc_box = GcFrame(self.open_sequences)
        rc_box = RcFrame(self.open_sequences)
        tr_box = TrFrame(self.open_sequences)
        
        nbook.append_page(os_box)
        nbook.append_page(en_box)
        nbook.append_page(sl_box)
        nbook.append_page(cm_box)
        nbook.append_page(gc_box)
        nbook.append_page(rc_box)
        nbook.append_page(tr_box)
        
        nbook.set_tab_label_text(os_box, "Open Sequences")
        nbook.set_tab_label_text(en_box, "Entrez Database")
        nbook.set_tab_label_text(sl_box, "Sequence Length")
        nbook.set_tab_label_text(cm_box, "Sequence Comparison")
        nbook.set_tab_label_text(gc_box, "GC Content")
        nbook.set_tab_label_text(rc_box, "Reverse Complement")
        nbook.set_tab_label_text(tr_box, "Translate Sequence")
        self.add(vbox)
        
        self.show_all()
        self.connect("delete-event", self.on_quit)
        

    def on_quit(self, widget, event):
        Gtk.main_quit()

       
window = MyWindow()

Gtk.main()
