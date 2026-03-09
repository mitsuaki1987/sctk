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
      $ git checkout qe-7.4.1
      $ git clone https://github.com/mitsuaki1987/sctk.git -b develop SCTK
      $ patch -p1 < SCTK/patch.diff

   開発版(developブランチ)を試す場合については https://github.com/mitsuaki1987/sctk/blob/develop/readme.md を参照
      
#. オリジナルの QuantumESPRESSO_ の場合と同様に ``configure`` で環境設定をする。
               
   .. code-block:: bash

       $ ./configure
       
#. メイク

   .. code-block:: bash

       $ make pw ph pp sctk
               
   実行ファイル名は ``sctk.x`` である。
   コンパイルが不完全に終わることがあり、その場合は再度 ``make`` を実行する。

   .. code-block:: bash

       $ make sctk
