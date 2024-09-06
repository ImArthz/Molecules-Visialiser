# Molecule Visualizer
## Introdução
O Molecule Visualizer é um aplicativo de visualização de moléculas que oferece uma interface gráfica para exibir moléculas em 2D e 3D. Utiliza a biblioteca RDKit para manipulação de moléculas e py3Dmol para visualizações 3D. Abaixo está uma explicação detalhada das funções e da estrutura do código.

## Executando o Código de Visualização de Moléculas
Este código é uma aplicação GUI (Interface Gráfica do Usuário) para visualizar moléculas usando Python. Abaixo, você encontrará uma explicação passo a passo sobre como executá-lo e as bibliotecas necessárias.

# Bibliotecas Necessárias
* RDKit: Biblioteca para manipulação e visualização de moléculas químicas.
* py3Dmol: Biblioteca para visualização 3D de moléculas.
* Tkinter: Biblioteca padrão do Python para criar interfaces gráficas.
* Pillow: Biblioteca para manipulação de imagens.
# Para instalar essas bibliotecas, use os seguintes comandos:

```bash
pip install rdkit
pip install py3Dmol
pip install pillow
```
* Nota: O RDKit pode ser difícil de instalar diretamente via pip, pois pode precisar de compilação. Em alguns casos, é melhor usar uma distribuição como o Anaconda, que já inclui o RDKit.

## Passos para Executar o Código
# Preparação do Ambiente

Crie um diretório de trabalho:

```bash
mkdir molecule_visualizer
cd molecule_visualizer
```
# Execute o Código

* Abra um terminal e navegue até o diretório onde você salvou o arquivo Python.

* Execute o script Python com o comando:

```bash
python main.py
```
A GUI (Interface Gráfica do Usuário) será exibida. Nela, você poderá buscar moléculas pelo nome e visualizá-las em 2D e 3D.
##Estrutura do Código
## 1. Importação de Módulos
O código começa com a importação dos módulos necessários:

```bash

import json
import os
import logging
import webbrowser
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import py3Dmol
import tkinter as tk
```
from tkinter import ttk, filedialog, messagebox
json: Manipula arquivos JSON para carregar e salvar dados.
os: Realiza operações do sistema operacional, como criação de diretórios e manipulação de arquivos.
logging: Registra informações sobre a execução do aplicativo, útil para depuração.
webbrowser: Abre URLs no navegador padrão.
rdkit: Biblioteca para química computacional; Chem e AllChem são usados para manipulação e visualização de moléculas.
py3Dmol: Para visualização 3D de moléculas.
tkinter: Biblioteca para criação da interface gráfica do usuário (GUI).
## 2. Função ensure_output_directory()
```bash

def ensure_output_directory():
    if not os.path.exists('outputs'):
        os.makedirs('outputs')
```
Objetivo: Garante que o diretório outputs exista.
Vantagem: Evita erros ao tentar salvar arquivos em um diretório que não existe.
## 3. Função load_molecules()
```bash

def load_molecules(file_path='database/molecules.json'):
    with open(file_path, 'r') as file:
        molecules = json.load(file)
    return molecules
```
Objetivo: Carrega as informações das moléculas de um arquivo JSON.
Formato Esperado do JSON:

```bash

{
    "Water": "O",
    "Methane": "C",
    "Ethanol": "CCO",
    "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
}
```
Vantagem: Permite que o aplicativo carregue dados de moléculas de forma estruturada e flexível.
## 4 . Função visualize_molecule()
```bash
def visualize_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logging.error("Invalid SMILES string: %s", smiles)
        return None
    AllChem.Compute2DCoords(mol)
    logging.info("Molecule visualized: %s", smiles)
    return mol
```
Objetivo: Cria a estrutura 2D da molécula a partir da string SMILES.
Vantagem: Gera uma representação visual básica da molécula, útil para a visualização rápida.
## 5. Função draw_molecule()
```bash
def draw_molecule(molecule, output_file='molecule.png'):
    ensure_output_directory()
    img = Draw.MolToImage(molecule)
    img.save(f'outputs/{output_file}')
    logging.info("Image saved as outputs/%s", output_file)
```
Objetivo: Desenha a molécula em 2D e salva como um arquivo PNG.
Vantagem: Permite a visualização em 2D da molécula, útil para análises visuais rápidas e relatórios.
## 6. Função visualize_with_py3dmol()
```bash

def visualize_with_py3dmol(smiles, output_file='molecule_3d.html'):
    ensure_output_directory()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logging.error("Invalid SMILES string: %s", smiles)
        return
    
    mol_3d = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_3d)
    AllChem.UFFOptimizeMolecule(mol_3d)
    
    mol_block = Chem.MolToMolBlock(mol_3d)
    
    viewer = py3Dmol.view(width=800, height=600)
    viewer.addModel(mol_block, "mol")
    viewer.setStyle({'stick': {}})
    
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

    html_file = f'outputs/{output_file}'
    with open(html_file, 'w') as f:
        f.write(viewer._make_html())
    logging.info("3D view saved as %s", html_file)
    return output_file
```
Objetivo: Gera uma visualização 3D da molécula e salva como um arquivo HTML.
Vantagem: Proporciona uma visualização interativa e detalhada da estrutura molecular, permitindo manipulação e exploração em 3D.
## 7. Função cleanup_files()
```bash
def cleanup_files(directory='outputs'):
    for filename in os.listdir(directory):
        if filename.endswith('.html') or filename.endswith('.png'):
            file_path = os.path.join(directory, filename)
            os.remove(file_path)
            logging.info("Deleted %s", file_path)
```
Objetivo: Remove arquivos antigos do diretório outputs.
Vantagem: Mantém o diretório limpo e evita acúmulo desnecessário de arquivos.
## 8. Função update_molecule_list()
```bash
def update_molecule_list():
    molecules = load_molecules()
    molecule_listbox.delete(0, tk.END)
    for name in molecules.keys():
        molecule_listbox.insert(tk.END, name)
```
Objetivo: Atualiza a lista de moléculas na GUI com base no arquivo JSON.
Vantagem: Mantém a interface gráfica atualizada com as moléculas disponíveis para visualização.
##9. Função on_search_change()
```bash
def on_search_change(*args):
    search_term = search_var.get().lower()
    filtered_molecules = [name for name in all_molecules if search_term in name.lower()]
    molecule_listbox.delete(0, tk.END)
    for name in filtered_molecules:
        molecule_listbox.insert(tk.END, name)
```
Objetivo: Filtra a lista de moléculas com base no termo de busca inserido pelo usuário.
Vantagem: Facilita a localização de moléculas específicas na lista, melhorando a usabilidade.
## 10. Função visualize_molecule_gui()
```bash
def visualize_molecule_gui():
    selected_molecule = molecule_listbox.get(tk.ACTIVE)
    if not selected_molecule:
        messagebox.showerror("Error", "Please select a molecule.")
        return

    molecules = load_molecules()
    if selected_molecule in molecules:
        smiles = molecules[selected_molecule]
        molecule_name = selected_molecule

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
```
Objetivo: Lida com a visualização de uma molécula selecionada na GUI. Gera imagens 2D e 3D da molécula e abre a visualização 3D no navegador.
Vantagem: Integra a funcionalidade de visualização com a interface gráfica, proporcionando uma experiência de usuário fluida.
## 11. Função open_git_link()
```bash
def open_git_link():
    webbrowser.open("https://github.com/creators")
```
Objetivo: Abre um link para o repositório Git dos criadores.
Vantagem: Permite aos usuários acessar rapidamente o repositório para mais informações ou contribuições.
## 12. Função exit_gui()
```bash
def exit_gui():
    cleanup_files()
    root.quit()
```
Objetivo: Limpa arquivos temporários e fecha a aplicação.
Vantagem: Garante que os arquivos temporários sejam removidos e que o aplicativo seja fechado corretamente.

## 13. Função main_gui()
```bash
def main_gui():
    global molecule_listbox, root, search_var, all_molecules
    root = tk.Tk()
    root.title("Molecule Visualizer")

    root.configure(bg='#e0f7fa')  

    tk.Label(root, text="Search for a molecule:", bg='#e0f7fa', font=('Helvetica', 14, 'bold')).pack(pady=10)

    search_var = tk.StringVar()
    search_var.trace_add("write", on_search_change)
    search_entry = tk.Entry(root, textvariable=search_var, width=80, font=('Helvetica', 12))
    search_entry.pack(pady=10)

    tk.Label(root, text="Select a molecule to visualize:", bg='#e0f7fa', font=('Helvetica', 14, 'bold')).pack(pady=10)

    molecule_listbox = tk.Listbox(root, width=80, height=20, bg='#ffffff', borderwidth=2, relief='solid', font=('Helvetica', 12))
    molecule_listbox.pack(pady=10, fill=tk.BOTH, expand=True)

    all_molecules = load_molecules().keys()
    update_molecule_list()

    button_frame = tk.Frame(root, bg='#e0f7fa')
    button_frame.pack(pady=10)

    tk.Button(button_frame, text="Visualize Molecule", command=visualize_molecule_gui, bg='#b3e5fc', fg='#000000', relief='raised', font=('Helvetica', 12)).pack(side=tk.LEFT, padx=5)
    tk.Button(button_frame, text="Exit", command=exit_gui, bg='#b3e5fc', fg='#000000', relief='raised', font=('Helvetica', 12)).pack(side=tk.LEFT, padx=5)
    tk.Button(button_frame, text="Git Link", command=open_git_link, bg='#b3e5fc', fg='#000000', relief='raised', font=('Helvetica', 12)).pack(side=tk.LEFT, padx=5)

    root.mainloop()
```
Objetivo: Cria e configura a interface gráfica principal. Define os widgets da GUI, como labels, entradas de texto, listas e botões.
Vantagem: Oferece uma interface intuitiva para interação com o aplicativo, facilitando a visualização de moléculas e a navegação.
## Conclusão
O Molecule Visualizer proporciona uma maneira poderosa e interativa de visualizar moléculas em 2D e 3D. Cada função foi projetada para modularizar e otimizar o processo de visualização e interação com as moléculas. Com uma interface gráfica amigável e recursos avançados de visualização, o aplicativo é uma ferramenta valiosa para pesquisadores e entusiastas da química.

<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.3/css/all.min.css">

## Dados

**Autor:** Arthur De Oliveira Mendonça 

**Redes Sociais:**

<div align="center">
    <h3>Connect with me:</h3>
    <p><strong>Arthur Mendonça</strong></p>
    <p>
        <a href="https://github.com/ImArthz" target="blank">
            <img align="center" src="https://raw.githubusercontent.com/rahuldkjain/github-profile-readme-generator/master/src/images/icons/Social/github.svg" alt="ImArthz" height="30" width="40" />
        </a>
        <a href="https://twitter.com/Im_Arthz" target="blank">
            <img align="center" src="https://raw.githubusercontent.com/rahuldkjain/github-profile-readme-generator/master/src/images/icons/Social/twitter.svg" alt="Im_Arthz" height="30" width="40" />
        </a>
        <a href="https://api.whatsapp.com/send?phone=37988528423" target="blank">
            <img align="center" src="https://raw.githubusercontent.com/rahuldkjain/github-profile-readme-generator/master/src/images/icons/Social/whatsapp.svg" alt="WhatsApp" height="30" width="40" />
        </a>
        <a href="https://discordapp.com/users/imarthz" target="blank">
            <img align="center" src="https://raw.githubusercontent.com/rahuldkjain/github-profile-readme-generator/master/src/images/icons/Social/discord.svg" alt="imarthz" height="30" width="40" />
        </a>
    </p>
</div>




