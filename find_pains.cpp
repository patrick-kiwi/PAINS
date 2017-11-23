#include <iostream>
#include <map>
#include <string>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <fstream>
#include <regex>

using namespace std;
using namespace OpenBabel;

map<string, string> get_smarts_from_file (string);
map<int, string> get_smiles_from_file (string);



int main(int argc,char **argv) {
    
    OBConversion conv;
    conv.SetInFormat("smiles");
    OBMol mol;
    
//Load the SMILES string of the parent molecules and 3 smarts patterns files into memory

    map<int, string> molmap = get_smiles_from_file ( argv[1] );
    map<string, string> ffa = get_smarts_from_file( argv[2] ); 
    map<string, string> ffb = get_smarts_from_file( argv[3] );
    map<string, string> ffc = get_smarts_from_file( argv[4] );



//iterate over target molecules and SMARTS fingerprints

    for (map<int, string>::iterator itt = molmap.begin(); itt != molmap.end(); itt++) { 
static int c = 0;
//cout << "Input smiles: " << itt->second << "\n";
        conv.ReadString(&mol, itt->second);
        
    
        for (map<string, string>::iterator it = ffa.begin(); it != ffa.end(); it++) {
            OBSmartsPattern sm;  //moving this outside the loop created a memory leak
            sm.Init(it->second);  //capture the smatrs match
            if  ( sm.IsValid() and sm.HasMatch(mol) ) 
                cout << "Found a match: " << it->first << "  Parent mol: " << itt->second << "\n";
        }
        
        
         for (map<string, string>::iterator it = ffb.begin(); it != ffb.end(); it++) {
            OBSmartsPattern sm;  //moving this outside the loop created a memory leak
            sm.Init(it->second);  //capture the smatrs match
            if  ( sm.IsValid() and sm.HasMatch(mol) ) 
                cout << "Found a match: " << it->first << "  Parent mol: " << itt->second << "\n";
        }
        
        
        
         for (map<string, string>::iterator it = ffc.begin(); it != ffc.end(); it++) {
            OBSmartsPattern sm;  //moving this outside the loop created a memory leak
            sm.Init(it->second);  //capture the smatrs match
            if  ( sm.IsValid() and sm.HasMatch(mol) ) 
                cout << "Found a match: " << it->first << "  Parent mol: " << itt->second << "\n";
        }
        
     c+=1;
     if (c>=5000)
         break;
    }


  return(0);
}

//function definitions

map<string, string> get_smarts_from_file (string input) {
    map<string, string> res;
    string eachLine;
    ifstream theFile ( input );
    const regex smarts_regex{"(FRAGMENT\\s+regId=(\\S+)\\(\\d+\\)\\s+(\\S+).*)"};
    if (!theFile.is_open() )
        cout << "Error, couldnt open smarts pattern file\n";
    while (getline (theFile, eachLine) ) {
        smatch match;
        regex_search(eachLine, match, smarts_regex);
        res[match[2]] = match[3];
    }
    theFile.close();
    return res;
}

map<int, string> get_smiles_from_file (string input) {
    static int count = 0;
    map<int, string> res;
    string eachLine;
    ifstream theFile ( input );
    const regex r{"(\\d+)\\s+(.*)"};
    if (!theFile.is_open() )
        cout << "Error, couldnt open smarts pattern file\n";
     while (getline (theFile, eachLine) ) {
        smatch match;
        regex_search(eachLine, match, r);
        res[count++] = match[2];
    }
    theFile.close();
    return res;
}
