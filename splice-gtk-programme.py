import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


class Splice(Gtk.Bin):
    def __init__(self):
        Gtk.Bin.__init__(self)

        self.builder = Gtk.Builder()
        self.builder.add_from_file("splice-page.glade")
        self.sp_box = self.builder.get_object("Splice-box")
        self.add(self.sp_box)

        self.entry = self.builder.get_object("Splice-entry")
        self.button = self.builder.get_object("Splice-button")
        self.label = self.builder.get_object("Splice-result")
        self.button.connect("clicked", self.clicked_callback)
        self.entry.connect("activate", self.enter_callback)

    def splice(self):
        entry_text = Seq(self.entry.get_text(), IUPAC.unambiguous_dna)
        records = list(SeqIO.parse(entry_text, "fasta"))
        #print (records)
        #ros_2295 = entry_text[0]
        #print(ros_2295.seq)
        print("blank line...")

        for rec in records:
            temp = str(ros_2295.seq).replace(str(rec.seq), "")
            ros_2295.seq = temp

        sp_result = (ros_2295.seq)
        self.label.set_text(str(sp_result))

    def clicked_callback(self, button):
        self.splice()

    def enter_callback(self, widget):
        self.splice()

class MyWindow(Gtk.Window):
    def __init__(self):
        Gtk.Window.__init__(self)
        self.set_size_request(300, 200)

        vbox = Gtk.VBox()
        self.set_title("RNA Splicing")

        sp_box = Splice()
        vbox.add(sp_box)
        self.add(vbox)

        self.show_all()
        self.connect("delete_event", self.on_quit)

    def on_quit(self, widget, event):
        Gtk.main_quit()

window = MyWindow()

Gtk.main()
