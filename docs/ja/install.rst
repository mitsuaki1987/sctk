.. _QuantumESPRESSO: https://www.quantum-espresso.org/resources/users-manual

インストール方法
================
 
要件
----

使用要件はオリジナルの QuantumESPRESSO_ と同じである。
          
インストール手順
----------------

#. 以下のようにして、QuantumESPRESSO_ を展開した中でSCTKのリポジトリをクローンする。

   .. code-block:: bash

      $ git clone https://gitlab.com/QEF/q-e.git
      $ cd q-e
      $ git checkout 96cdd5ac6af9c060be392a95f14dbcbca5c1a890
      $ git clone https://github.com/mitsuaki1987/sctk.git -b sctk1.2.1-qe6.7
      $ patch -p1 < sctk/patch.diff

   開発版を試す場合には次のようにする。

   .. code-block:: bash

      $ git clone https://gitlab.com/QEF/q-e.git
      $ cd q-e
      $ git checkout 0e7fd34c024bde2ef2545e3aeda2bbfec770d500
      $ git clone https://github.com/mitsuaki1987/sctk.git -b develop
      $ patch -p1 < sctk/patch.diff
      
#. オリジナルの QuantumESPRESSO_ の場合と同様に ``configure`` で環境設定をする。
               
   .. code-block:: bash

       $ ./configure
       
#. メイク

   .. code-block:: bash

       $ make pw ph pp
       $ cd sctk
       $ make
               
   実行ファイル名は ``sctk.x`` である。

