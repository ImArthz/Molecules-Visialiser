import json
import os
import logging
import webbrowser
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import py3Dmol
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

# Configuração do logging
logging.basicConfig(filename='app.log', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Garantir que o diretório de saída exista
def ensure_output_directory():
    if not os.path.exists('outputs'):
        os.makedirs('outputs')

# Carregar moléculas do arquivo JSON
def load_molecules(file_path='database/molecules.json'):
    with open(file_path, 'r') as file:
        molecules = json.load(file)
    return molecules

# Gerar estrutura 2D da molécula
def visualize_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logging.error("Invalid SMILES string: %s", smiles)
        return None
    AllChem.Compute2DCoords(mol)
    logging.info("Molecule visualized: %s", smiles)
    return mol

# Desenhar molécula e salvar como PNG
def draw_molecule(molecule, output_file='molecule.png'):
    ensure_output_directory()
    img = Draw.MolToImage(molecule)
    img.save(f'outputs/{output_file}')
    logging.info("Image saved as outputs/%s", output_file)

# Visualizar molécula em 3D e salvar como HTML
def visualize_with_py3dmol(smiles, output_file='molecule_3d.html'):
    ensure_output_directory()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logging.error("Invalid SMILES string: %s", smiles)
        return
    
    # Gerar coordenadas 3D
    mol_3d = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_3d)
    AllChem.UFFOptimizeMolecule(mol_3d)
    
    # Converter molécula para coordenadas 3D
    mol_block = Chem.MolToMolBlock(mol_3d)
    
    # Usar py3Dmol para visualizar a molécula em 3D com rótulos para átomos
    viewer = py3Dmol.view(width=800, height=600)
    viewer.addModel(mol_block, "mol")
    viewer.setStyle({'stick': {}})
    
    # Adicionar rótulos para cada átomo na molécula
    for atom in mol_3d.GetAtoms():
        pos = mol_3d.GetConformer().GetAtomPosition(atom.GetIdx())
        viewer.addLabel(atom.GetSymbol(), {
            'position': {'x': pos.x, 'y': pos.y, 'z': pos.z},
            'backgroundColor': 'white',
            'fontColor': 'black',
            'fontSize': 14,
            'showBackground': True
        })
    
    viewer.zoomTo()

    # Salvar visualização 3D em um arquivo HTML
    html_file = f'outputs/{output_file}'
    with open(html_file, 'w') as f:
        f.write(viewer._make_html())
    logging.info("3D view saved as %s", html_file)
    return html_file

# Gerar e salvar a fórmula química como imagem PNG
def draw_formula(molecule, output_file='formula.png'):
    ensure_output_directory()
    img = Draw.MolToImage(molecule, size=(300, 300))
    img.save(f'outputs/{output_file}')
    logging.info("Formula image saved as outputs/%s", output_file)

# Limpar todos os arquivos HTML e PNG
def cleanup_files(directory='outputs'):
    for filename in os.listdir(directory):
        if filename.endswith('.html') or filename.endswith('.png'):
            file_path = os.path.join(directory, filename)
            os.remove(file_path)
            logging.info("Deleted %s", file_path)

# Atualizar lista de moléculas na GUI
def update_molecule_list():
    molecules = load_molecules()
    molecule_listbox.delete(0, tk.END)
    for name in molecules.keys():
        molecule_listbox.insert(tk.END, name)

# Função de busca com auto-completamento
def on_search_change(*args):
    search_term = search_var.get().lower()
    filtered_molecules = [name for name in all_molecules if search_term in name.lower()]
    molecule_listbox.delete(0, tk.END)
    for name in filtered_molecules:
        molecule_listbox.insert(tk.END, name)

# Funções para a GUI
def visualize_molecule_gui():
    selected_molecule = molecule_listbox.get(tk.ACTIVE)
    if not selected_molecule:
        messagebox.showerror("Error", "Please select a molecule.")
        return

    molecules = load_molecules()
    if selected_molecule in molecules:
        smiles = molecules[selected_molecule]
        molecule_name = selected_molecule

        # Limpar arquivos antigos
        cleanup_files()

        mol = visualize_molecule(smiles)
        if mol:
            draw_molecule(mol, f"{molecule_name}.png")
            html_file = visualize_with_py3dmol(smiles, f"{molecule_name}_3d.html")

            if html_file:
                file_path = os.path.join('outputs', html_file)
                webbrowser.open(file_path)
            else:
                messagebox.showerror("Error", "Failed to create HTML file.")
    else:
        messagebox.showerror("Error", "Selected molecule not found in the database.")

def visualize_formula_gui():
    selected_molecule = molecule_listbox.get(tk.ACTIVE)
    if not selected_molecule:
        messagebox.showerror("Error", "Please select a molecule.")
        return

    molecules = load_molecules()
    if selected_molecule in molecules:
        smiles = molecules[selected_molecule]
        molecule_name = selected_molecule

        # Limpar arquivos antigos
        cleanup_files()

        mol = visualize_molecule(smiles)
        if mol:
            draw_formula(mol, f"{molecule_name}_formula.png")
            img_file = os.path.join('outputs', f"{molecule_name}_formula.png")
            # Abrir imagem PNG gerada
            webbrowser.open(img_file)
        else:
            messagebox.showerror("Error", "Failed to visualize molecule.")

def search_molecule_online():
    selected_molecule = molecule_listbox.get(tk.ACTIVE)
    if not selected_molecule:
        messagebox.showerror("Error", "Please select a molecule.")
        return

    molecules = load_molecules()
    if selected_molecule in molecules:
        smiles = molecules[selected_molecule]
        # URL do PubChem para pesquisa de moléculas
        pubchem_url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{selected_molecule.replace(' ', '_')}"
        webbrowser.open(pubchem_url)
    else:
        messagebox.showerror("Error", "Selected molecule not found in the database.")

def open_git_link():
    webbrowser.open("https://github.com/ImArthz/Molecules-Visualiser")

def exit_gui():
    # Limpar arquivos e fechar o aplicativo
    cleanup_files()
    root.quit()

def main_gui():
    global molecule_listbox, root, search_var, all_molecules
    root = tk.Tk()
    root.title("Molecule Visualizer")

    # Configurações de estilo
    root.configure(bg='#e0f7fa')  # Fundo azul claro

    tk.Label(root, text="Search for a molecule:", bg='#e0f7fa', font=('Helvetica', 14, 'bold')).pack(pady=10)

    search_var = tk.StringVar()
    search_var.trace_add("write", on_search_change)
    search_entry = tk.Entry(root, textvariable=search_var, width=80, font=('Helvetica', 12))
    search_entry.pack(pady=10)

    tk.Label(root, text="Select a molecule to visualize:", bg='#e0f7fa', font=('Helvetica', 14, 'bold')).pack(pady=10)

    # Frame para a Listbox com a barra de rolagem
    listbox_frame = tk.Frame(root)
    listbox_frame.pack(pady=10)

    scrollbar = tk.Scrollbar(listbox_frame)
    scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

    molecule_listbox = tk.Listbox(listbox_frame, width=80, height=20, bg='#ffffff', borderwidth=2, relief='solid', font=('Helvetica', 12), yscrollcommand=scrollbar.set)
    molecule_listbox.pack(side=tk.LEFT, fill=tk.BOTH)

    scrollbar.config(command=molecule_listbox.yview)

    # Carregar moléculas e preencher a lista
    all_molecules = load_molecules().keys()
    update_molecule_list()

    button_frame = tk.Frame(root, bg='#e0f7fa')
    button_frame.pack(pady=10)

    tk.Button(button_frame, text="Visualize Molecule", command=visualize_molecule_gui, bg='#b3e5fc', fg='#000000', relief='raised', font=('Helvetica', 12)).pack(side=tk.LEFT, padx=5)
    tk.Button(button_frame, text="Visualize 2D", command=visualize_formula_gui, bg='#b3e5fc', fg='#000000', relief='raised', font=('Helvetica', 12)).pack(side=tk.LEFT, padx=5)
    tk.Button(button_frame, text="Search Online", command=search_molecule_online, bg='#b3e5fc', fg='#000000', relief='raised', font=('Helvetica', 12)).pack(side=tk.LEFT, padx=5)
    tk.Button(button_frame, text="Exit", command=exit_gui, bg='#b3e5fc', fg='#000000', relief='raised', font=('Helvetica', 12)).pack(side=tk.LEFT, padx=5)
    tk.Button(button_frame, text="Git Link", command=open_git_link, bg='#b3e5fc', fg='#000000', relief='raised', font=('Helvetica', 12)).pack(side=tk.LEFT, padx=5)

    root.mainloop()

if __name__ == "__main__":
    main_gui()
