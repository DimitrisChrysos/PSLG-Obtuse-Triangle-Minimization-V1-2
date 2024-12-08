- Περιγραφή του προγράμματος:
    - Πιστεύουμε πως το πρόγραμμα μας υλοποιεί όλες τις απαιτήσεις της εργασίας που μας ζητήθηκαν μέσω της εκφώνησης, των φροντιστηρίων και των συζητήσεων στο eclass.
    - Το ex.cpp είναι το αρχείο με την main που καλεί τις απαραίτητες συναρτήσεις για την υλοποίηση της εργασίας μας.
    - Αρχικά διαβάζουμε κάποιο json αρχείο ως input
    - Κατασκευάζουμε το cdt με τις οδηγίες που μας δίνονται
    - Διαλέγουμε τη μέθοδο (ή τις μεθόδους) που θα χρησιμοποιηθούν και τις παραμέτρους τους, από το αρχείο που διαβάστηκε
    - Βάζουμε τα απαραίτητα steiner points για να μηδενίσουμε όσο γίνεται τα αμβλυγώνια τρίγωνα.
    - Δημιουργούμε το output αρχείο json.



- Κατάλογος των αρχείων κώδικα και το ρόλο τους:
1) Αρχείο ex.cpp στο βασικό directory:
    1) Εκτελεί την παραπάνω συνοπτική περιγραφή που αναφέρθηκε.
    2) Χρησιμοποιεί της συναρτήσεις:
        - handle_methods(), ανάλογα τις παραμέτρους επιλέγει τη μέθοδο και το τροόπο που θα χρησιμοποιηθεί:
            1) local_search()
            2) sim_annealing()
            3) ant_colony_optimization()
            - Κάθε μέθοδος χρησιμοποιεί βοηθητικές συναρτήσεις και δομές
        - write_output(), δημιουργεί το output json αρχείο
2) Κατάλογος includes:
    1) Υποκατάλογος custom_cdt_class:
        1. Αρχείο custom_cdt_class.hpp, περιέχει την custom κλάση που δημιουργήθηκε για τις ανάγκες των μεθόδων:
            - insert_no_flip()
            - insert_steiner_x_y(), αποθηκεύει σε ένα vertex στο cdt τις συντεταγμένες των steiner points που εισάγουμε
            - my_is_flippable(), συμβάλει στον έλεγχο για το αν μία ακμή έχει τη δυνατότητα να γίνει flipped
            - remove_no_flip()
            - Ορίζει κάποια typedef για τις ανάγκες όλων των αρχείων της εργασίας
        2. Αρχείο custom_cdt_class.cpp, βοηθητικό για την μεταγλώττιση
    2) Υποκατάλογος read_write_file:
        1. Αρχείο read_write_file.hpp, περιέχει τη δήλωση όλων των βοηθητικών συναρτήσεων για την ανάγνωση του αρχείου json και για τη δημιουργία του output  αρχείο json
        2. Αρχείο read_write_file.cpp, περιέχει την υλοποίηση όλων των βοηθητικών συναρτήσεων για την ανάγνωση του αρχείου json και για τη δημιουργία του output  αρχείο json
    3) Υποκατάλογος steiner_methods:
        1. Αρχείο steiner_methods.hpp περιέχει τη δήλωση των 5 μεθόδων steiner και μίας βοηθητικής
        2. Αρχείο steiner_methods.cpp περιέχει την υλοποίηση των 5 μεθόδων steiner και μίας βοηθητικής
    4) Υποκατάλογος utils:
        1. Αρχείο utils.hpp περιέχει τη δήλωση των περισσότερων βοηθητικών μεθόδων και κλάσεων
        2. Αρχείο utils.hpp περιέχει την υλοποίηση των περισσότερων βοηθητικών μεθόδων και κλάσεων
3) Κατάλογος outputs:
    - Περιέχει και για τις τρεις μεθόδους (ls, sa, ant) τα outputs των test που μας δόθηκαν
4) Κατάλογοι test:
    - Περιέχουν τα test που μας δόθηκαν για τις τρεις μεθόδους



- Οδηγίες μεταγλώττισης:
    1) cmake -DCGAL_DIR=/usr/lib/CGAL
    2) make



- Οδηγίες χρήσης του προγράμματος:
    1) Για το example τρέχουμε:
        - ls: ./opt_triangulation -i input_ls.json -o outputs/output_ls.json
        - sa: ./opt_triangulation -i input_sa.json -o outputs/output_sa.json
        - ant: ./opt_triangulation -i input_ant.json -o outputs/output_ant.json
    2) Για τα tests:
        - ls: ./opt_triangulation -i tests_ls/instance_1.json -o outputs/tests_ls/output1.json
        - sa: ./opt_triangulation -i tests_sa/instance_1.json -o outputs/tests_sa/output1.json
        - ant: ./opt_triangulation -i tests_ant/instance_1.json -o outputs/tests_ant/output1.json
        - Αντίστοιχα αλλάζουμε του αριθμούς 1, για όποιο test θέλουμε
    3) Για τα test_instances:
        - ls: ./opt_triangulation -i test_instances_ls/instance_test_1.json -o outputs/test_instances_ls/output1.json
        - sa: ./opt_triangulation -i test_instances_sa/instance_test_1.json -o outputs/test_instances_sa/output1.json
        - ant: ./opt_triangulation -i test_instances_ant/instance_test_1.json -o outputs/test_instances_ant/output1.json
        - Αντίστοιχα αλλάζουμε του αριθμούς 1, για όποιο test θέλουμε



- Πλήρη στοιχεία των φοιτητών που ανέπτυξαν την εργασία:
1) Ονοματεπώνυμο: Αναστάσιος Μουμουλίδης - ΑΜ: 1115202100108
2) Ονοματεπώνυμο: Δημήτριος Χρυσός - ΑΜ: 1115202100275
    


- Επιλογή παραμέτρων:
    1) Η επιλογή των παραμέτρων, έγινε για κάθε input json αρχείο και μέθοδο ξεχωριστά
    2) Η επιλογή έγινε με βάση δικές μας δοκιμές και τα αποτελέσματα που βγάζαμε



- Σύγκριση μεθόδων:
    1) Με βάση τα αρχεία output που αναφερθήκαμε προηγουμένως, είναι φανερό ότι η πιο αποτελεσματική μέθοδος είναι η local search, μετά η simulated annealing και μετά η ant colony optimization
    2) Επίσης παρατηρείται ότι η χρήση της παραμέτρου, "delaunay": false, βοηθάει μερικές φορές κυρίως τις μεθόδους sa και ants αλλά χωρίς σταθερά αποτελέσματα για κάθε input



- github repo link: https://github.com/DimitrisChrysos/project-emiris.git



