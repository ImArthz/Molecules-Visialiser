# Molecule Visualizer

## Introdução

O Molecule Visualizer é um aplicativo de visualização de moléculas que oferece uma interface gráfica para exibir moléculas em 2D e 3D. Utiliza a biblioteca RDKit para manipulação de moléculas e py3Dmol para visualizações 3D. Este README fornece uma visão geral das funcionalidades e da estrutura do código.

## Bibliotecas Necessárias

### Para executar o Molecule Visualizer, você precisará das seguintes bibliotecas:

- **RDKit**: Biblioteca para manipulação e visualização de moléculas químicas.
- **py3Dmol**: Biblioteca para visualização 3D de moléculas.
- **Tkinter**: Biblioteca padrão do Python para criar interfaces gráficas.
- **Pillow**: Biblioteca para manipulação de imagens.

### Para instalar essas bibliotecas, use os seguintes comandos:

```bash
pip install rdkit
pip install py3Dmol
pip install pillow
```
* Nota: O RDKit pode ser difícil de instalar diretamente via pip, pois pode precisar de compilação. Em alguns casos, é melhor usar uma distribuição como o Anaconda, que já inclui o RDKit.

## Passos para Executar o Código
* Preparação do Ambiente
* Crie um diretório de trabalho:
```bash
mkdir molecule_visualizer
cd molecule_visualizer
```
* Execute o Código
* Abra um terminal e navegue até o diretório onde você salvou o arquivo Python.
* Execute o script Python com o comando:
```bash
python main.py
```
* A GUI (Interface Gráfica do Usuário) será exibida. Nela, você poderá buscar moléculas pelo nome e visualizá-las em 2D e 3D.

##  Estrutura do Código
# 1. Importação de Módulos
```bash

import json
import os
import logging
import webbrowser
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import py3Dmol
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
```
## Funções do Código
* Abaixo estão descrições detalhadas das principais funções do código:

### ensure_output_directory()
```bash
# Garantir que o diretório de saída exista
def ensure_output_directory():
    if not os.path.exists('outputs'):
        os.makedirs('outputs')
```

* Descrição: Garante que o diretório outputs exista para salvar as saídas do programa.
* Funcionamento: Verifica se o diretório 'outputs' existe. Caso contrário, ele o cria.

### load_molecules(file_path='database/molecules.json')
```bash
# Carregar moléculas do arquivo JSON
def load_molecules(file_path='database/molecules.json'):
    with open(file_path, 'r') as file:
        molecules = json.load(file)
    return molecules
```

* Descrição: Carrega as informações das moléculas de um arquivo JSON.
* Funcionamento: Lê o arquivo JSON especificado e retorna um dicionário com as moléculas e suas respectivas strings SMILES.
  
### visualize_molecule(smiles)
```bash
# Gerar estrutura 2D da molécula
def visualize_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logging.error("Invalid SMILES string: %s", smiles)
        return None
    AllChem.Compute2DCoords(mol)
    logging.info("Molecule visualized: %s", smiles)
    return mol
```

* Descrição: Gera a estrutura 2D da molécula a partir de uma string SMILES.
* Funcionamento: Usa o RDKit para converter a string SMILES em um objeto de molécula e calcula suas coordenadas 2D.
  
### draw_molecule(molecule, output_file='molecule.png')
```bash
# Desenhar molécula e salvar como PNG
def draw_molecule(molecule, output_file='molecule.png'):
    ensure_output_directory()
    img = Draw.MolToImage(molecule)
    img.save(f'outputs/{output_file}')
    logging.info("Image saved as outputs/%s", output_file)
```

* Descrição: Desenha a molécula em 2D e salva como um arquivo PNG.
* Funcionamento: Utiliza a função MolToImage do RDKit para gerar a imagem da molécula e salva no diretório 'outputs'.

### visualize_with_py3dmol(smiles, output_file='molecule_3d.html')
```bash
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
```

* Descrição: Gera uma visualização 3D da molécula e salva como um arquivo HTML.
* Funcionamento: Converte a molécula para um bloco de coordenadas 3D, utiliza o py3Dmol para criar a visualização e salva como HTML.
  
### draw_formula(molecule, output_file='formula.png')
```bash
# Gerar e salvar a fórmula química como imagem PNG
def draw_formula(molecule, output_file='formula.png'):
    ensure_output_directory()
    img = Draw.MolToImage(molecule, size=(300, 300))
    img.save(f'outputs/{output_file}')
    logging.info("Formula image saved as outputs/%s", output_file)
```
* Descrição: Gera e salva a fórmula química da molécula como imagem PNG.
* Funcionamento: Usa o RDKit para desenhar a fórmula da molécula e salva como um arquivo PNG no diretório 'outputs'.
  
### cleanup_files(directory='outputs')
```bash
# Limpar todos os arquivos HTML e PNG
def cleanup_files(directory='outputs'):
    for filename in os.listdir(directory):
        if filename.endswith('.html') or filename.endswith('.png'):
            file_path = os.path.join(directory, filename)
            os.remove(file_path)
            logging.info("Deleted %s", file_path)
```
* Descrição: Limpa arquivos antigos do diretório outputs.
* Funcionamento: Remove todos os arquivos HTML e PNG existentes no diretório especificado.
  
### update_molecule_list()
```bash
# Atualizar lista de moléculas na GUI
def update_molecule_list():
    molecules = load_molecules()
    molecule_listbox.delete(0, tk.END)
    for name in molecules.keys():
        molecule_listbox.insert(tk.END, name)
```
* Descrição: Atualiza a lista de moléculas na GUI com base no arquivo JSON.
* Funcionamento: Recarrega o arquivo de moléculas e atualiza a Listbox da interface com os nomes das moléculas disponíveis.
### on_search_change(*args)
```bash
# Função de busca com auto-completamento
def on_search_change(*args):
    search_term = search_var.get().lower()
    filtered_molecules = [name for name in all_molecules if search_term in name.lower()]
    molecule_listbox.delete(0, tk.END)
    for name in filtered_molecules:
        molecule_listbox.insert(tk.END, name)
```
* Descrição: Filtra a lista de moléculas com base no termo de busca inserido pelo usuário.
* Funcionamento: Monitora mudanças no campo de busca e atualiza a lista de moléculas exibidas na GUI conforme o usuário digita.
  
### open_in_browser(url)
```bash
# Função para abrir URL em navegadores especificados
def open_in_browser(url):
    # Lista de navegadores a tentar
    browsers = {
        'chrome': 'google-chrome',
        'firefox': 'firefox',
        'safari': 'safari',
        'opera': 'opera'
    }
    
    # Tenta abrir o URL com os navegadores especificados
    for name, command in browsers.items():
        try:
            webbrowser.get(command).open(url)
            return
        except webbrowser.Error:
            continue

    # Se todos falharem, abre no navegador padrão
    webbrowser.open(url)
```
* Descrição: Abre uma URL no navegador web.
* Funcionamento: Tenta abrir a URL com navegadores específicos (Chrome, Firefox, Safari, Opera). Caso falhe, abre no navegador padrão do sistema.
  
### visualize_molecule_gui()
```bash
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
                file_path = os.path.join('outputs', f"{molecule_name}_3d.html")
                open_in_browser(file_path)  # Abre o arquivo HTML com a função
            else:
                messagebox.showerror("Error", "Failed to create HTML file.")
    else:
        messagebox.showerror("Error", "Selected molecule not found in the database.")
```
* Descrição: Lida com a visualização de uma molécula selecionada na GUI. Gera imagens 2D e 3D da molécula e abre a visualização 3D no navegador.
* Funcionamento: Busca a molécula selecionada pelo usuário, gera e exibe as visualizações apropriadas, e exibe mensagens de erro caso algo falhe.
  
### visualize_formula_gui()
```bash
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
            open_in_browser(img_file)
        else:
            messagebox.showerror("Error", "Failed to visualize molecule.")
```
* Descrição: Gera e exibe a fórmula 2D da molécula selecionada na GUI.
* Funcionamento: Similar à visualize_molecule_gui, mas foca na fórmula da molécula.

### search_molecule_online()
```bash
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
        open_in_browser(pubchem_url)
    else:
        messagebox.showerror("Error", "Selected molecule not found in the database.")
```
* Descrição: Busca informações detalhadas sobre a molécula online.
* Funcionamento: Abre uma página do PubChem com informações detalhadas da molécula selecionada.

### open_git_link()
```bash
def open_git_link():
    open_in_browser("https://github.com/ImArthz/Molecules-Visualiser")
```
* Descrição: Abre o repositório GitHub do projeto.
* Funcionamento: Usa open_in_browser para abrir o link do repositório.

### exit_gui()
```bash
def exit_gui():
    # Limpar arquivos e fechar o aplicativo
    cleanup_files()
    root.quit()
```
* Descrição: Limpa arquivos temporários e fecha a aplicação.
* Funcionamento: Chama cleanup_files() para limpar a pasta 'outputs' e encerra a GUI.
  
### main_gui()
```bash
def main_gui():
    global molecule_listbox, root, search_var, all_molecules
    ensure_output_directory()  # Garante que a pasta de saída exista
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
```
* Descrição: Cria e configura a interface gráfica principal.
* Funcionamento: Configura a janela principal da aplicação, define os componentes da interface como campos de busca, lista de moléculas, e botões de interação.

## Novidades
* Adição do Botão "Chemistry Resource"
* O aplicativo agora inclui um botão "Chemistry Resource" que abre um site renomado para uma explicação detalhada da molécula pesquisada.

## Como usar:
* Digite o nome da molécula no campo de busca.
* Selecione a molécula da lista exibida.
* Clique no botão "Chemistry Resource" para ser redirecionado a uma página com informações detalhadas sobre a molécula.
* Nota: O link para a página de química detalhada é um exemplo e pode precisar ser ajustado conforme a disponibilidade dos recursos.

## Conclusão
O Molecule Visualizer proporciona uma maneira poderosa e interativa de visualizar moléculas em 2D e 3D. Cada função foi projetada para modularizar e otimizar o processo de visualização e interação com as moléculas. Com uma interface gráfica amigável e recursos avançados de visualização, o aplicativo é uma ferramenta valiosa para pesquisadores e entusiastas da química.

# Colaboradores
* Arthur De Oliveira Mendonça

### Redes Sociais:

<div align="center"> <h3>Connect with me:</h3> <p><strong>Arthur Mendonça</strong></p> <p> <a href="https://github.com/ImArthz" target="blank"> <img align="center" src="https://raw.githubusercontent.com/rahuldkjain/github-profile-readme-generator/master/src/images/icons/Social/github.svg" alt="ImArthz" height="30" width="40" /> </a> <a href="https://twitter.com/Im_Arthz" target="blank"> <img align="center" src="https://raw.githubusercontent.com/rahuldkjain/github-profile-readme-generator/master/src/images/icons/Social/twitter.svg" alt="Im_Arthz" height="30" width="40" /> </a> <a href="https://api.whatsapp.com/send?phone=37988528423" target="blank"> <img align="center" src="https://raw.githubusercontent.com/rahuldkjain/github-profile-readme-generator/master/src/images/icons/Social/whatsapp.svg" alt="WhatsApp" height="30" width="40" /> </a> <a href="https://discordapp.com/users/imarthz" target="blank"> <img align="center" src="https://raw.githubusercontent.com/rahuldkjain/github-profile-readme-generator/master/src/images/icons/Social/discord.svg" alt="imarthz" height="30" width="40" /> </a> </p> </div>
